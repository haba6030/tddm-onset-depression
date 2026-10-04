#!/usr/bin/env bash
# Regenerate every reported statistic and figure from the shipped fits. No model fitting.
# About 20 minutes; the starting-point simulation takes most of it.
set -euo pipefail
cd "$(dirname "$0")"
[ -f data/study4_HC.csv ] || python data/fetch_data.py
CANON=results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv
mkdir -p results/logs
step () {                    # step <label> <command...>; output goes to results/logs/<label>.log
  local log="results/logs/$(echo "$1" | tr -c 'A-Za-z0-9\n-' '_').log"; shift
  echo "== $log"; "$@" > "$log" 2>&1 || { echo "FAILED, see $log"; exit 1; }
}

step "C++ vs R likelihood"     Rscript scripts/tDDM/verify_cpp_likelihood.R
step "group statistics"        Rscript scripts/tDDM/group_stats.R --fits "$CANON"
step "main-model reanalysis"   Rscript scripts/tDDM/memo_0919/reanalysis_accurate.R
step "model-free behaviour"    Rscript scripts/tDDM/memo_0919/model_free_bestrun.R
step "demographics (Table 1)"  Rscript scripts/tDDM/demographics_table.R --fits "$CANON"
step "CES-D continuous"        Rscript scripts/tDDM/ces_continuous_correlation.R
step "S1 collinearity"         Rscript scripts/paper/s1_collinearity.R
step "S2 optimizer runs"       Rscript scripts/paper/s2_optimizer_runs_accurate.R
step "controls summary"        Rscript scripts/tDDM/memo_0919/stage6_summary_accurate.R
step "supplement tables"       Rscript scripts/paper/supplement_tables_accurate.R
step "RT signature"            Rscript scripts/tDDM/rt_signature_model_check.R --fits "$CANON" --cores 6
step "starting-point export"   Rscript scripts/tDDM/memo_0919/starting_point_export.R
step "starting-point bound"    python scripts/tDDM/memo_0919/starting_point_sensitivity.py

step "figures 1-2"  python scripts/visualization/make_schematic_figures.py
step "figures 3-5"  python scripts/visualization/make_paper_figures_v2.py
step "figures 6-7"  python scripts/visualization/make_boundary_figures.py
step "figure S1"    python scripts/tDDM/figS1_individual_recovery.py
step "figure S2"    python scripts/tDDM/figS2_absolute_fit.py
step "figures S3-4" python scripts/eeg/figs_cesd_and_age.py
echo "done: logs in results/logs/, figures in docs/PAPER/figures/"
