#!/usr/bin/env bash
# Reproduce the paper. Run from anywhere; paths are relative to the repository root.
#
#   bash reproduce.sh regenerate                 statistics and figures from the shipped fits (~20 min)
#   bash reproduce.sh fits       [--cores N]     refit the baseline, contrast and time-varying models (~2 h)
#   bash reproduce.sh validation [--cores N] [--fits FILE]
#                                                identifiability and validation analyses (~14 h)
#   bash reproduce.sh checks                     checks of the likelihood and the simulator (minutes)
#
# Options: --cores N (default 8), --fits FILE (validation only; default: the paper's fits CANON),
#          --dry-run (print the commands without running them).
# Environment: RSCRIPT and PYTHON choose the interpreters (default Rscript, python).
# Each call writes results/logs/<stage>_<time>/ with one log per step and run_info.txt
# (command, commit, versions, start and end time, exit code).
set -euo pipefail
cd "$(dirname "$0")"

usage () { sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

STAGE="${1:-}"; [ -n "$STAGE" ] || usage 1; shift
CORES=8; DRY=0
CANON=results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv
FITS="$CANON"
while [ $# -gt 0 ]; do
  case "$1" in
    --cores) CORES="$2"; shift 2 ;;
    --fits) FITS="$2"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) usage 0 ;;
    *) echo "unknown option: $1"; usage 1 ;;
  esac
done
case "$STAGE" in regenerate|fits|validation|checks) ;; -h|--help) usage 0 ;; *) echo "unknown stage: $STAGE"; usage 1 ;; esac
R="${RSCRIPT:-Rscript}"; PY="${PYTHON:-python}"

if [ "$DRY" = 1 ]; then
  step () { shift; echo "$*"; }
else
  [ -f data/study4_HC.csv ] && [ -f data/study4_dep.csv ] || "$PY" data/fetch_data.py
  LOGDIR="results/logs/${STAGE}_$(date +%Y%m%d_%H%M%S)"
  [ -e "$LOGDIR" ] && { echo "$LOGDIR exists; wait a second and rerun"; exit 1; }
  mkdir -p "$LOGDIR"
  {
    echo "command: bash reproduce.sh $STAGE --cores $CORES --fits $FITS"
    if git rev-parse --short HEAD > /dev/null 2>&1; then
      echo "commit: $(git rev-parse --short HEAD)$(git diff HEAD --quiet || echo ' (uncommitted changes)')"
    else echo "commit: unknown (not a git checkout)"; fi
    echo "R: $("$R" --version 2>/dev/null | head -1)"
    echo "python: $("$PY" --version 2>&1)"
    echo "start: $(date '+%Y-%m-%d %H:%M:%S')"
  } > "$LOGDIR/run_info.txt"
  trap 'code=$?; { echo "end: $(date "+%Y-%m-%d %H:%M:%S")"; echo "exit_code: $code"; } >> "$LOGDIR/run_info.txt"' EXIT
  n=0
  step () {                    # step <label> <command...>; output goes to $LOGDIR/<nn>_<label>.log
    n=$((n + 1)); local log; log="$LOGDIR/$(printf '%02d' "$n")_$1.log"; shift
    echo "== $(date '+%H:%M:%S') $*"; echo "$*" >> "$LOGDIR/commands.txt"
    "$@" > "$log" 2>&1 || { echo "FAILED, see $log"; exit 1; }
  }
fi

latest () { ls -dt results/tDDM/"$1"_* | head -1; }

case "$STAGE" in
regenerate)
  step cpp_likelihood        "$R" scripts/tDDM/verify_cpp_likelihood.R
  step model_free            "$R" scripts/tDDM/memo_0919/model_free_bestrun.R
  step baseline_model        "$R" scripts/paper/baseline_primary_accurate.R
  step group_stats           "$R" scripts/tDDM/group_stats.R --fits "$CANON"
  step reanalysis            "$R" scripts/tDDM/memo_0919/reanalysis_accurate.R
  step demographics          "$R" scripts/tDDM/demographics_table.R --fits "$CANON"
  step cesd_continuous       "$R" scripts/tDDM/ces_continuous_correlation.R
  step s1_collinearity       "$R" scripts/paper/s1_collinearity.R
  step s2_optimizer_runs     "$R" scripts/paper/s2_optimizer_runs_accurate.R
  step controls_summary      "$R" scripts/tDDM/memo_0919/stage6_summary_accurate.R
  step supplement_tables     "$R" scripts/paper/supplement_tables_accurate.R
  step rt_signature          "$R" scripts/tDDM/rt_signature_model_check.R --fits "$CANON" --cores 6
  step starting_point_export "$R" scripts/tDDM/memo_0919/starting_point_export.R
  step starting_point_bound  "$PY" scripts/tDDM/memo_0919/starting_point_sensitivity.py
  export FIG_GROUP_LABELS='Low CES-D,High CES-D'        # group names printed in the paper figures
  step fig1        "$PY" scripts/visualization/make_schematic_figures.py --merged
  step fig2_4      "$PY" scripts/visualization/make_paper_figures_v2.py
  step figS2_S3    "$PY" scripts/visualization/make_boundary_figures.py
  step figS1       "$PY" scripts/tDDM/figS1_individual_recovery.py
  step figS4       "$PY" scripts/tDDM/figS2_absolute_fit.py
  step figS5_S6    "$PY" scripts/eeg/figs_cesd_and_age.py
  ;;
fits)
  step baseline_model  "$R" scripts/tDDM/refit_corrected_likelihood.R --arms 4a_static1
  step contrast_model  "$R" scripts/tDDM/refit_corrected_likelihood.R --arms 4c_static
  # time-varying model: 10 runs per signal ordering, then the best run per participant
  step tv_mean_first      "$R" scripts/tDDM/optimizer_reproducibility.R --reps 10 --cores "$CORES"
  step tv_mean_first_best "$R" scripts/tDDM/memo_0919/build_bestof_runs.R \
    --runs-dir "$( [ "$DRY" = 1 ] && echo '<latest optimizer_reproducibility_mean_r1_*>' || latest optimizer_reproducibility_mean_r1)" \
    --out results/tDDM/refit_bestof10/fits_4j_mean_first_best10.csv
  step tv_contrast_first      "$R" scripts/tDDM/optimizer_reproducibility.R --orientation contrast --reps 10 --cores "$CORES"
  step tv_contrast_first_best "$R" scripts/tDDM/memo_0919/build_bestof_runs.R \
    --runs-dir "$( [ "$DRY" = 1 ] && echo '<latest optimizer_reproducibility_contrast_r1_*>' || latest optimizer_reproducibility_contrast_r1)" \
    --out results/tDDM/refit_bestof10/fits_4j_contrast_first_best10.csv
  ;;
validation)
  # Generating onsets for the controls: the sample mean and the low depressive-symptom group mean
  # of the fitted onsets, with tau <= 0 set to 0 (0.891 and 0.804 for the paper's fits).
  read -r TAU_COMMON TAU_BASE < <(FITS_FILE="$FITS" "$R" -e 'f <- read.csv(Sys.getenv("FITS_FILE"));
    t <- pmax(f$tau, 0); cat(sprintf("%.3f %.3f\n", mean(t), mean(t[f$group == "HC"])))' | tail -1) || true
  [ -n "${TAU_BASE:-}" ] || { echo "could not read onsets from $FITS"; exit 1; }
  echo "fits: $FITS (tau-common $TAU_COMMON, tau-base $TAU_BASE)"
  DE="--de-iter 200 --de-np 80 --restarts 3 --cores $CORES"
  step recovery_tv        "$R" scripts/tDDM/recovery_redesign_tau.R --trans image --arms redesign --trials 20,80 \
    --N 50 --reps 3 --restarts 3 --cores "$CORES"
  step recovery_baseline  "$R" scripts/tDDM/recovery_baseline.R --N 50 --reps 3 --restarts 3 --trials 20,80 --cores "$CORES"
  step split_half_mean    "$R" scripts/tDDM/split_half_predictive.R --K 30 --models 4a,4j_mean_first --trans image \
    --itermax 100 --np 60 --restarts 2 --cores "$CORES"
  step split_half_contrast "$R" scripts/tDDM/split_half_predictive.R --K 30 --models 4j_contrast_first --trans image \
    --itermax 100 --np 60 --restarts 2 --cores "$CORES"
  step split_half_threshold "$R" scripts/tDDM/split_half_predictive.R --K 30 --models 4a --cores "$CORES"
  step negative_control   "$R" scripts/tDDM/nulltau_corrected.R --reps 8 $DE --fits "$FITS" --tau-common "$TAU_COMMON"
  step positive_control_20 "$R" scripts/tDDM/positive_control_tau.R --reps 10 --gaps 0,0.085,0.169,0.34 --trans image \
    --tau-base "$TAU_BASE" $DE --fits "$FITS"
  step positive_control_80 "$R" scripts/tDDM/positive_control_tau.R --reps 5 --gaps 0,0.34 --ntrials 80 --trans image \
    --tau-base "$TAU_BASE" $DE --fits "$FITS"
  step rt_signature       "$R" scripts/tDDM/rt_signature_model_check.R --fits "$FITS" --cores "$CORES"
  step choice_signature   "$R" scripts/tDDM/memo_0919/choice_signature_check.R --cores "$CORES"
  step starting_point_005 "$R" scripts/tDDM/nulltau_corrected.R --reps 5 $DE --fits "$FITS" --tau-common "$TAU_COMMON" --wgap 0.05
  step starting_point_010 "$R" scripts/tDDM/nulltau_corrected.R --reps 5 $DE --fits "$FITS" --tau-common "$TAU_COMMON" --wgap 0.10
  step ndt_bounds         "$R" scripts/tDDM/nDT_bound_ladder_engine.R --bounds 1.5,2.0,2.5,3.0,3.5 $DE
  step tau_causes_rt      "$R" scripts/tDDM/tau_causes_rt_check.R --clamp-first --cores "$CORES"
  ;;
checks)
  step transition_density   "$R" scripts/tDDM/verify_transition_density.R
  step cpp_likelihood       "$R" scripts/tDDM/verify_cpp_likelihood.R
  step simulator            "$R" scripts/tDDM/diag_simulator_vs_likelihood.R
  step quadrature_spectral  "$R" scripts/tDDM/diag_quadrature_spectral.R
  ;;
esac

[ "$DRY" = 1 ] || echo "done: logs in $LOGDIR/"
