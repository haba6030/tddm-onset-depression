#!/usr/bin/env bash
# 2026-09-28: re-run every control under the accurate split-at-tau quadrature (engine
# default TDDM_QUAD = "accurate") and the bridge-corrected simulator. Settings match
# the grid-era stage-1/2 runs (run_stage1*.sh, run_stage2.sh) except the operating
# points, which come from the accurate canonical fits (tau <= 0 set to 0):
#   pooled mean tau 0.891 (null-tau common onset), HC mean tau 0.804 (positive-control base).
# Order = manuscript priority (main-text items first).
set -u
cd "$(dirname "$0")/../../.."
M=results/tDDM/stage6_accurate
FITS=${FITS:-results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv}
R=${RESTARTS:-3}
C=${CORES:-12}
mkdir -p "$M"
run () {                      # run <name> <command...>; skips steps already done
  local name=$1; shift
  if [ -f "$M/$name.done" ]; then echo "[skip ] $name" >> "$M/runner.log"; return; fi
  echo "[start] $name  $(date '+%F %T')" >> "$M/runner.log"
  if "$@" > "$M/$name.log" 2>&1; then
    date '+%F %T' > "$M/$name.done"; echo "[done ] $name  $(date '+%F %T')" >> "$M/runner.log"
  else
    echo "[FAIL ] $name  $(date '+%F %T') -- see $M/$name.log" >> "$M/runner.log"
  fi
}
# --- simulation only (no refits) ---
run choicesig Rscript scripts/tDDM/memo_0919/choice_signature_check.R --cores "$C"
run rtsig     Rscript scripts/tDDM/rt_signature_model_check.R --fits "$FITS" --cores "$C"
run taurt     Rscript scripts/tDDM/tau_causes_rt_check.R --cores "$C"
# --- refits ---
run contrastfirst Rscript scripts/tDDM/optimizer_reproducibility.R --orientation contrast \
  --reps 10 --cores "$C"
run splithalf Rscript scripts/tDDM/split_half_predictive.R \
  --K 30 --models 4a,4j_mean_first --trans image \
  --itermax 100 --np 60 --restarts 2 --cores "$C"
run splithalf_contrast Rscript scripts/tDDM/split_half_predictive.R \
  --K 30 --models 4j_contrast_first --trans image \
  --itermax 100 --np 60 --restarts 2 --cores "$C"
run nulltau Rscript scripts/tDDM/nulltau_corrected.R \
  --reps 8 --restarts "$R" --de-iter 200 --de-np 80 --cores "$C" \
  --fits "$FITS" --tau-common 0.891
run recovery Rscript scripts/tDDM/recovery_redesign_tau.R --trans image --arms redesign \
  --trials 20,80 --N 50 --reps 3 --restarts "$R" --cores "$C"
run ndtladder Rscript scripts/tDDM/nDT_bound_ladder_engine.R \
  --bounds 1.5,2.0,2.5,3.0,3.5 --de-iter 200 --de-np 80 --restarts "$R" --cores "$C"
for gap in 0.05 0.10; do
  run "wgap$gap" Rscript scripts/tDDM/nulltau_corrected.R \
    --reps 5 --restarts "$R" --de-iter 200 --de-np 80 --cores "$C" \
    --fits "$FITS" --tau-common 0.891 --wgap "$gap"
done
run poscontrol Rscript scripts/tDDM/positive_control_tau.R \
  --reps 10 --gaps 0,0.085,0.169,0.34 --trans image --tau-base 0.804 \
  --de-iter 200 --de-np 80 --restarts "$R" --cores "$C" --fits "$FITS"
run ladder80 Rscript scripts/tDDM/positive_control_tau.R \
  --reps 5 --gaps 0,0.34 --ntrials 80 --trans image --tau-base 0.804 \
  --de-iter 200 --de-np 80 --restarts "$R" --cores "$C" --fits "$FITS"
echo "[stage 6 finished] $(date '+%F %T')" >> "$M/runner.log"
