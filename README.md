# Depressive symptoms, evidence accumulation and caution in a social decision: code and results

This repository reproduces the analyses in *Depressive symptoms are associated with weaker evidence accumulation and greater caution during social decision-making* (Kim & Yao). The paper fits drift-diffusion models to ultimatum-game choices and response times from Experiment 4 of Heffner et al. (2021). The sample is 351 adults: 205 scored below the CES-D cutoff of 16 (low depressive-symptom group) and 146 scored at or above it (high depressive-symptom group).

The primary analysis uses the time-stationary baseline model. In this model a mean reward–emotion signal drives evidence accumulation for the whole trial, and its drift rate and decision threshold are compared between the groups. The paper then compares the baseline model with two larger ones. The contrast model adds a weight for the reward–emotion contrast signal. The time-varying model lets the contrast signal enter only after an onset τ. BIC favours the baseline model. The group difference in τ is an exploratory analysis, reported together with the identifiability and control analyses that limit how it can be read.

The repository contains:

- the likelihood engine (`models/tDDM/tddm_engine.R`, with a C++ version of the time-varying likelihood in `models/tDDM/tddm_split.cpp`) and every script that produces a reported statistic, table, or figure;
- per-participant parameter estimates for all models and all optimizer runs;
- the outputs of the identifiability and validation analyses (Supplement S4) and of the non-decision-time bound analysis (Supplement S3);
- checks of the likelihood and of the trial simulator;
- a script that downloads the behavioural data from the original authors and rebuilds the exact files we analysed.

## Names used in the code

The scripts and output files use the internal names on the right. Group outputs are reported as `DEP − HC`, which is the high minus the low depressive-symptom group. `HC` and `DEP` are legacy labels for the two CES-D groups, which are defined by a self-report scale.

| Paper | Code and output files |
|---|---|
| baseline model (mean-signal drift, threshold, non-decision time) | `4a`, `4a_static1` |
| contrast model (adds the contrast weight) | `4c`, `4c_static` |
| time-varying model (adds the onset τ; mean signal first) | `4j`, `4j_mean_first` |
| time-varying model with the contrast signal first | `4j_contrast_first`, `--orientation contrast` |
| low depressive-symptom group (CES-D < 16) | `HC` |
| high depressive-symptom group (CES-D ≥ 16) | `DEP` |
| mean-signal drift, contrast weight, threshold, non-decision time | `d_mean`, `d_contrast`, `thres`, `nDT` |

## Requirements

- R ≥ 4.3 with `DEoptim`, `Rcpp`, `lme4`, `BayesFactor`, and a C++ compiler that `Rcpp::sourceCpp` can use (Xcode command-line tools on macOS, `build-essential` on Linux). Fitting uses `parallel::mclapply`, which forks, so the fitting scripts run on macOS or Linux and not on Windows.
- Python ≥ 3.10 with `numpy`, `pandas`, `scipy`, `matplotlib`, `statsmodels`.

We ran the fits and most analyses with R 4.3.3 (DEoptim 2.2.8, Rcpp 1.1.0, lme4 1.1.37, BayesFactor 0.9.12-4.7) and Python 3.12 (numpy 1.26.4, pandas 2.2.2, scipy 1.17.1, matplotlib 3.10.7, statsmodels 0.14.6). The three runs added in November 2026 (baseline recovery, baseline threshold split-half, onset shift and response time) used R 4.4.1 with DEoptim 2.2.8. Each of them records its command and versions in `run_info.txt`.

Run every command from the repository root. The first R script that sources the engine compiles `tddm_split.cpp` into `models/tDDM/.rcpp_cache/`, which takes about a minute once.

## 1. Data

The raw data are not redistributed here because the source repository carries no licence. Download and rebuild them:

```bash
python data/fetch_data.py
```

The script downloads `study4_ug.csv` and `study4_ec.csv` from [jpheffne/epe](https://github.com/jpheffne/epe) at a pinned commit. It then applies the original study's preregistered exclusion: a participant is kept when their rating of "neutral" falls inside the 100 × 100 px square at the centre of the emotion grid, which leaves 351 of 508. Each signal is divided by its SD over the retained trials of all participants. The SD is the population SD, and the signals are not centred. The script writes `data/study4_HC.csv` and `data/study4_dep.csv` and checks both against the SHA-256 hashes of the files the paper was fitted on. It stops if either hash differs.

The model reads `RPE_scale` and `EPE_X1_scale` (valence prediction error). Response times outside 0.15–10 s are dropped inside the engine (`tddm_load_study4`), which removes 281 of 7,020 trials.

## 2. Regenerate statistics and figures from the shipped fits (about 20 minutes)

These commands use only the files in `results/` and run no model fits. Figures go to `docs/PAPER/figures/`. Section 4 lists what each command reproduces.

```bash
CANON=results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv
Rscript scripts/tDDM/verify_cpp_likelihood.R
Rscript scripts/tDDM/memo_0919/model_free_bestrun.R
Rscript scripts/paper/baseline_primary_accurate.R
Rscript scripts/tDDM/group_stats.R --fits "$CANON"
Rscript scripts/tDDM/memo_0919/reanalysis_accurate.R
Rscript scripts/tDDM/demographics_table.R --fits "$CANON"
Rscript scripts/tDDM/ces_continuous_correlation.R
Rscript scripts/paper/s1_collinearity.R
Rscript scripts/paper/s2_optimizer_runs_accurate.R
Rscript scripts/tDDM/memo_0919/stage6_summary_accurate.R
Rscript scripts/paper/supplement_tables_accurate.R
Rscript scripts/tDDM/rt_signature_model_check.R --fits "$CANON" --cores 6
Rscript scripts/tDDM/memo_0919/starting_point_export.R
python scripts/tDDM/memo_0919/starting_point_sensitivity.py

export FIG_GROUP_LABELS='Low CES-D,High CES-D'                 # group names printed in the paper figures
python scripts/visualization/make_schematic_figures.py --merged   # Fig. 1
python scripts/visualization/make_paper_figures_v2.py             # Figs. 2-4
python scripts/visualization/make_boundary_figures.py             # Figs. S2-S3
python scripts/tDDM/figS1_individual_recovery.py                  # Fig. S1
python scripts/tDDM/figS2_absolute_fit.py                         # Fig. S4
python scripts/eeg/figs_cesd_and_age.py                           # Figs. S5-S6
```

The figure file names come from earlier drafts and do not follow the paper's numbering. The figure table in Section 4 gives the mapping.

## 3. Refit from scratch (about 15 hours on 8 cores)

Each fitting command writes a new timestamped directory under `results/tDDM/`. Set `C` to the number of cores available.

```bash
C=8
latest () { ls -dt results/tDDM/"$1"_* | head -1; }

# baseline and contrast models (<1 h each)
Rscript scripts/tDDM/refit_corrected_likelihood.R --arms 4a_static1
Rscript scripts/tDDM/refit_corrected_likelihood.R --arms 4c_static

# time-varying model: 10 runs per signal ordering, best run per participant (~20 min each)
Rscript scripts/tDDM/optimizer_reproducibility.R --reps 10 --cores $C
Rscript scripts/tDDM/memo_0919/build_bestof_runs.R --runs-dir "$(latest optimizer_reproducibility_mean_r1)" \
  --out results/tDDM/refit_bestof10/fits_4j_mean_first_best10.csv
Rscript scripts/tDDM/optimizer_reproducibility.R --orientation contrast --reps 10 --cores $C
Rscript scripts/tDDM/memo_0919/build_bestof_runs.R --runs-dir "$(latest optimizer_reproducibility_contrast_r1)" \
  --out results/tDDM/refit_bestof10/fits_4j_contrast_first_best10.csv
FITS=results/tDDM/refit_bestof10/fits_4j_mean_first_best10.csv

# identifiability and validation analyses (~14 h). --tau-common 0.891 and --tau-base 0.804 are the
# sample mean and the low depressive-symptom group mean of the paper's onsets with tau <= 0 set to 0;
# recompute them for new fits.
Rscript scripts/tDDM/recovery_redesign_tau.R --trans image --arms redesign --trials 20,80 --N 50 \
  --reps 3 --restarts 3 --cores $C
Rscript scripts/tDDM/recovery_baseline.R --N 50 --reps 3 --restarts 3 --trials 20,80 --cores $C
Rscript scripts/tDDM/split_half_predictive.R --K 30 --models 4a,4j_mean_first --trans image \
  --itermax 100 --np 60 --restarts 2 --cores $C
Rscript scripts/tDDM/split_half_predictive.R --K 30 --models 4j_contrast_first --trans image \
  --itermax 100 --np 60 --restarts 2 --cores $C
Rscript scripts/tDDM/split_half_predictive.R --K 30 --models 4a --cores $C
Rscript scripts/tDDM/nulltau_corrected.R --reps 8 --restarts 3 --de-iter 200 --de-np 80 --cores $C \
  --fits "$FITS" --tau-common 0.891
Rscript scripts/tDDM/positive_control_tau.R --reps 10 --gaps 0,0.085,0.169,0.34 --trans image \
  --tau-base 0.804 --de-iter 200 --de-np 80 --restarts 3 --cores $C --fits "$FITS"
Rscript scripts/tDDM/positive_control_tau.R --reps 5 --gaps 0,0.34 --ntrials 80 --trans image \
  --tau-base 0.804 --de-iter 200 --de-np 80 --restarts 3 --cores $C --fits "$FITS"
Rscript scripts/tDDM/rt_signature_model_check.R --fits "$FITS" --cores $C
Rscript scripts/tDDM/memo_0919/choice_signature_check.R --cores $C
for gap in 0.05 0.10; do
  Rscript scripts/tDDM/nulltau_corrected.R --reps 5 --restarts 3 --de-iter 200 --de-np 80 --cores $C \
    --fits "$FITS" --tau-common 0.891 --wgap $gap
done
Rscript scripts/tDDM/nDT_bound_ladder_engine.R --bounds 1.5,2.0,2.5,3.0,3.5 --de-iter 200 --de-np 80 \
  --restarts 3 --cores $C
Rscript scripts/tDDM/tau_causes_rt_check.R --clamp-first --cores $C

# checks of the likelihood and the simulator (minutes each)
Rscript scripts/tDDM/verify_transition_density.R
Rscript scripts/tDDM/verify_cpp_likelihood.R
Rscript scripts/tDDM/diag_simulator_vs_likelihood.R
Rscript scripts/tDDM/diag_quadrature_spectral.R
```

The summary and figure scripts in Section 2 read the shipped run directories. To use new fits, point them at the new directories (Section 4).

Refits do not give bit-identical parameter estimates. The fitting scripts use `RNGkind("L'Ecuyer-CMRG")` with `mclapply`, so the random stream each participant receives depends on how participants are distributed over cores. DEoptim also resolves near-ties between candidate solutions differently when two likelihood implementations differ in the tenth decimal place. The group-level results are stable. Supplement S2 reports the spread across the ten runs of the time-varying model, and Section 5 compares the shipped fits with a refit through the public fitting path.

## 4. Where each result comes from

The Results report, in order: model-free behaviour, the model comparison, the baseline model's drift and threshold, and the exploratory onset analysis with its controls.

`BASE` = `results/tDDM/corrected_likelihood_20260819_012638/fits_4a_static1.csv`, the baseline model's per-participant estimates. It has no onset, so its likelihood does not involve the split-at-τ integral described in Section 6.

`CANON` = `results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv`, the time-varying model's per-participant estimates (mean signal first). For each participant it keeps the run with the lowest negative log-likelihood among 10 optimizer runs. The column `tau` holds the fitted value. Every analysis sets τ ≤ 0 to 0, because any τ ≤ 0 gives the same single-phase model.

Runtimes are wall-clock times on 8 to 12 cores (Apple M-series).

### Model fits

| Paper item | Script and command | Output | Time |
|---|---|---|---|
| Baseline model = `BASE` | `Rscript scripts/tDDM/refit_corrected_likelihood.R --arms 4a_static1` | `BASE` | <1 h |
| Contrast model | `Rscript scripts/tDDM/refit_corrected_likelihood.R --arms 4c_static` | `results/tDDM/corrected_likelihood_20260819_005659/fits_4c_static.csv` | <1 h |
| 10 optimizer runs of the time-varying model, seeds 20260821–30 (the runs behind `CANON`) | `Rscript scripts/tDDM/diag_quadrature_refit.R --seed <seed> --cores 12`, once per seed | `results/tDDM/diag_quadrature_refit_20260927_*/fits_4j_accurate_rep01seed.csv` | ~26 min per run |
| `CANON` (best of 10) | `Rscript scripts/tDDM/memo_0919/build_accurate_bestrun.R` | `CANON`, `summary.txt` | s |
| 10 runs with the contrast signal first | `Rscript scripts/tDDM/optimizer_reproducibility.R --orientation contrast --reps 10 --cores 8` | `results/tDDM/optimizer_reproducibility_contrast_r1_20260928_005730/` | ~20 min |
| Contrast-first best of 10 | `Rscript scripts/tDDM/memo_0919/build_bestof_runs.R --runs-dir <that folder> --out <file>` | `.../canonical_4j_accurate_bestof10/fits_4j_accurate_contrast_first_best10.csv` | s |

### Main text, Table 1 and Supplement S1–S3, S5–S6

Paths starting with `.../` are inside `results/tDDM/`.

| Paper item | Script | Output |
|---|---|---|
| Model-free comparison: acceptance, response time, trial-level choice model | `Rscript scripts/tDDM/memo_0919/model_free_bestrun.R` | stdout |
| Model comparison: BIC ladder baseline → contrast → time-varying (Fig. 2) | `Rscript scripts/tDDM/memo_0919/reanalysis_accurate.R` (section 6) | `.../canonical_4j_accurate_bestof10/reanalysis_accurate.txt` |
| Model comparison: lowest-BIC model by group; BIC of the two signal orderings | `Rscript scripts/paper/baseline_primary_accurate.R` | `.../canonical_4j_accurate_bestof10/baseline_primary.txt` |
| Baseline model: group comparisons on drift and threshold (Fig. 3a–b); baseline rows of Tables S9 and S10 | `Rscript scripts/paper/baseline_primary_accurate.R` | `.../canonical_4j_accurate_bestof10/baseline_primary.txt` |
| Time-varying model: group comparisons on drift, contrast weight, threshold, non-decision time and τ with τ ≤ 0 set to 0 (Fig. 3c); bound rates; τ ≤ 0 share and robustness; CES-D, age, gender; covariate regression and its Bayes factor | `Rscript scripts/tDDM/memo_0919/reanalysis_accurate.R` (sections 1–5) | `.../canonical_4j_accurate_bestof10/reanalysis_accurate.txt` |
| Group comparisons for any fits file, on the values as fitted (τ is not set to 0, so its group means differ from the paper's) | `Rscript scripts/tDDM/group_stats.R --fits <file>` | stdout |
| Table 1 | `Rscript scripts/tDDM/demographics_table.R --fits CANON` | `results/tDDM/cohort_demographics/` |
| CES-D as a continuous measure | `Rscript scripts/tDDM/ces_continuous_correlation.R` | stdout |
| S1: collinearity of the two signals | `Rscript scripts/paper/s1_collinearity.R` | `results/paper_numbers/` |
| Table S2: transition-density check | `Rscript scripts/tDDM/verify_transition_density.R` | `.../transition_density_check_20260820_003321/` |
| S2: agreement among the 10 runs | `Rscript scripts/paper/s2_optimizer_runs_accurate.R` | stdout |
| Table S3: non-decision-time bounds 1.5–3.5 s | `nDT_bound_ladder_engine.R` (Section 3); summary by `supplement_tables_accurate.R` | `.../nDT_bound_ladder_engine_20260928_192521/` |
| Tables S9 and S10, time-varying rows | `Rscript scripts/paper/supplement_tables_accurate.R` | `.../canonical_4j_accurate_bestof10/supplement_tables.txt` |

### Identifiability and validation analyses (Supplement S4, Table S4)

The analysis numbers follow Table S4. All runs used the commands in Section 3. Each simulated participant was fitted three times and the lowest-NLL fit was kept, except in Analysis 2 (Supplement S2 gives its reduced budget). `supplement_tables_accurate.R` and `stage6_summary_accurate.R` summarize these runs into `supplement_tables.txt` and `stage6_summary.txt` in `canonical_4j_accurate_bestof10/`.

| Analysis | Paper item | Script | Output | Time |
|---|---|---|---|---|
| 1. Individual recovery, baseline model | Table S5a | `recovery_baseline.R` | `recovery_baseline_20261106_210914/` | 0.5 h |
| 1. Individual recovery, time-varying model | Tables S5b–d, Fig. S1 | `recovery_redesign_tau.R` | `recovery_redesign_tau_20260928_192519/` | 0.4 h |
| 2. Split-half stability: baseline drift, mean-first onset | Fig. 4 | `split_half_predictive.R --models 4a,4j_mean_first` | `split_half_predictive_image_20260928_115546/` | 1.9 h |
| 2. Split-half stability: contrast-first onset | Fig. 4 | `split_half_predictive.R --models 4j_contrast_first` | `split_half_predictive_image_20260928_134940/` | 0.7 h |
| 2. Split-half stability: baseline threshold | Fig. 4 | `split_half_predictive.R --models 4a` | `split_half_predictive_image_20261106_205634/` | 0.2 h |
| 3. Negative control, equal group onsets | Table S6, Fig. S2b | `nulltau_corrected.R --tau-common 0.891` | `nulltau_corrected_image_20260928_211818/` | 1.1 h |
| 4. Positive control, 20 trials, 4 levels | Table S7, Fig. S3a, latent onset difference (Fig. S3c) | `positive_control_tau.R --gaps 0,0.085,0.169,0.34` | `positive_control_tau_20260928_231637/` | 2.2 h |
| 4. Positive control, 80 trials | Fig. S3b | `positive_control_tau.R --gaps 0,0.34 --ntrials 80` | `positive_control_tau_20260929_012719/` | 2.9 h |
| 5. Absolute fit: response-time signature | Table S8, Fig. S4 | `rt_signature_model_check.R --fits CANON` | `rt_signature_check_20260928_005349/` | minutes |
| 5. Absolute fit: choice signature | S4 text | `memo_0919/choice_signature_check.R` | `choice_signature_check_20260928_005021/` | minutes |
| 6. Starting-point control, separations 0.05 and 0.10 | Table S4 | `nulltau_corrected.R --wgap 0.05` (and `0.10`) | `nulltau_corrected_image_20260928_222147/`, `..._224927/` | 0.5 h each |
| Starting-point bound from acceptance rates | S4 text | `memo_0919/starting_point_export.R`, then `memo_0919/starting_point_sensitivity.py` | `memo_analyses_0919/starting_point_sensitivity_<date>/` (written when run) | ~12 min |
| Separability from general slowing: covariate regression | Fig. S2a | `memo_0919/reanalysis_accurate.R` (section 5) | `canonical_4j_accurate_bestof10/reanalysis_accurate.txt` | s |
| Separability from general slowing: onset shift and median response time | S4 text | `tau_causes_rt_check.R --clamp-first` | `tau_causes_rt_clampfirst_20261107_104449/` | minutes |
| Contrast-weight subgroups | S4 text | `supplement_tables_accurate.R` | `supplement_tables.txt` | s |

Some block titles inside `supplement_tables.txt` follow an earlier numbering of the Supplement. Its "Table S4: individual recovery" block is Table S5b–d in the paper. The other block titles (Tables S3 and S6–S10) match the paper.

`--tau-common 0.891` and `--tau-base 0.804` are the sample mean and the low depressive-symptom group mean of the onsets in `CANON`, with τ ≤ 0 set to 0.

### Figures

All figure scripts write to `docs/PAPER/figures/`. Set `FIG_GROUP_LABELS='Low CES-D,High CES-D'` to print the group names used in the paper. `FIG_OUT=<dir>` writes to another folder.

| Paper figure | Script | File |
|---|---|---|
| Fig. 1 (schematic, no data) | `python scripts/visualization/make_schematic_figures.py --merged` | `F1_model_space` |
| Fig. 2 | `python scripts/visualization/make_paper_figures_v2.py` | `F4_model_comparison` |
| Fig. 3 | the same | `F3_group_differences` |
| Fig. 4 | the same | `F5_split_half` |
| Fig. S1 | `python scripts/tDDM/figS1_individual_recovery.py` | `figS1_individual_recovery` |
| Figs. S2, S3 | `python scripts/visualization/make_boundary_figures.py` | `F6_onset_controls`, `F7_onset_sensitivity` |
| Fig. S4 | `python scripts/tDDM/figS2_absolute_fit.py` (run `rt_signature_model_check.R` first; its trial-level file contains observed response times and is not shipped) | `figS2_absolute_fit` |
| Figs. S5, S6 | `python scripts/eeg/figs_cesd_and_age.py` | `figS3_cesd_continuous`, `figS4_age_control` |

## 5. Reproducing the time-varying fits through the public fitting path

`CANON` was fitted with `scripts/tDDM/diag_quadrature_refit.R`, which evaluates the time-varying likelihood in R. The refit commands in Section 3 fit through the engine instead (`optimizer_reproducibility.R`), which evaluates the same likelihood in C++ and is 10 to 17 times faster. The two implementations agree to within 1.3 × 10⁻¹⁰ in log-likelihood on the stored fits, on 3,000 random parameter points, and on simulated data (`verify_cpp_likelihood.R`, output in `canonical_4j_accurate_bestof10/verify_cpp_likelihood.txt`).

We refitted the time-varying model 10 times through the engine (`results/tDDM/optimizer_reproducibility_mean_r1_20261002_130237/`, 17 minutes on 12 cores) and took the best run per participant with `build_bestof_runs.R`:

| | Δτ (high − low) | p | mean ΔBIC (time-varying − contrast) | ΔBIC < 0 | ΔBIC < −2 |
|---|---|---|---|---|---|
| `CANON` (paper) | +0.209 s | .044 | +1.790 | 10.0% | 2.8% |
| engine refit | +0.226 s | .031 | +1.789 | 10.0% | 2.8% |

The best negative log-likelihood agreed to within 10⁻⁶ for 84.9% of participants and to within 0.01 for all but two, where the engine run found a lower value (by up to 0.125). Fitted onsets correlated at r = 0.983. The engine refit's Δτ lies inside the range the paper reports across its own ten runs (+0.171 to +0.280 s, Supplement S2).

## 6. Notes on the likelihood

- `models/tDDM/tddm_engine.R` is the only implementation of the likelihood. The time-varying likelihood integrates the state density at τ (method of images) times the second-phase first-passage density over the evidence level. The engine does this by adaptive Gauss–Kronrod quadrature (relative and absolute tolerance 10⁻⁶) with breakpoints at a distance of min(bsep/4, 20√(t − τ)) from each boundary, because the second-phase density concentrates near the boundaries when a response falls shortly after τ. Supplement S1 describes the method, and `diag_quadrature_spectral.R` checks it against an independent spectral solution.
- `TDDM_QUAD` in the engine also accepts `"grid"`, a fixed 30-point rule we used earlier. It overstates the likelihood of a response that falls just after τ, and no reported result uses it. It is kept so that the comparison in `diag_quadrature_refit_*/summary.txt` can be reproduced.
- The engine's `--trans pipeline` option keeps an earlier transition density that ignored absorption and mis-scaled the diffusion. It exists so that the density comparison in Table S2 can be reproduced. Every reported fit uses `--trans image`.
- The trial simulator used by the identifiability analyses (`tddm_simulate_trial`, `tddm_simulate_block`) steps in 1 ms, ends a path with the Brownian-bridge probability that it crossed a boundary within the step, and cuts the step at τ. `diag_simulator_vs_likelihood.R` checks its decision-time and choice distributions against the likelihood.
- DEoptim runs with its default value-to-reach, so every fit runs all 200 generations. Setting `VTR` in `DEoptim.control()` makes DEoptim stop at generation 0.

## Repository layout

| Folder | Contents |
|---|---|
| `models/tDDM/` | likelihood engine and its C++ part |
| `scripts/` | fitting, analysis and figure scripts ([`scripts/README.md`](scripts/README.md)) |
| `results/tDDM/` | per-participant fits and the run folders cited above ([`results/README.md`](results/README.md)) |
| `data/` | `fetch_data.py`; the downloaded CSV files are not tracked |
| `docs/PAPER/figures/` | output folder of the figure scripts; figures are not tracked |
| `tests/` | `python -m pytest tests/` checks that every path this README cites exists and that every script and run folder is documented |

## Citation

Please cite the paper and the original dataset:

Heffner, J., Son, J.-Y., & FeldmanHall, O. (2021). Emotion prediction errors guide socially adaptive behaviour. *Nature Human Behaviour*, 5, 1391–1401.

## Licence

Code: MIT (see `LICENSE`). The derived per-participant estimates and simulation outputs in `results/` are released under CC BY 4.0. The behavioural data belong to the original authors and are not included.
