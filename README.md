# Depression and the onset of information in a social decision: code and results

This repository reproduces the analyses in *Depression relates to when, and not only how strongly, information enters a social decision*. The paper fits a time-varying drift-diffusion model to ultimatum-game choices and response times from Experiment 4 of Heffner et al. (2021). In that model, a mean reward–emotion signal drives accumulation from the start of the trial, and a reward–emotion contrast signal enters after an onset τ. The model is fitted to 351 adults (205 below and 146 at or above the CES-D cutoff of 16).

The repository contains:

- the likelihood engine (`models/tDDM/tddm_engine.R`, with a C++ version of the time-varying likelihood in `models/tDDM/tddm_split.cpp`) and every script that produces a reported statistic, table, or figure;
- per-participant parameter estimates for all models and all optimizer runs;
- the outputs of the identifiability analyses (individual recovery, split-half replication, positive control, negative control, starting-point sensitivity, non-decision-time bounds, response-time and choice signatures);
- checks of the likelihood and of the trial simulator;
- a script that downloads the behavioural data from the original authors and rebuilds the exact files we analysed.

## Requirements

- R ≥ 4.3 with `DEoptim`, `Rcpp`, `lme4`, `BayesFactor`, and a C++ compiler that `Rcpp::sourceCpp` can use (Xcode command-line tools on macOS, `build-essential` on Linux). Fitting uses `parallel::mclapply`, which forks, so the fitting scripts run on macOS or Linux and not on Windows.
- Python ≥ 3.10 with `numpy`, `pandas`, `scipy`, `matplotlib`, `statsmodels`.

We ran everything with R 4.3.3 (DEoptim 2.2.8, Rcpp 1.1.0, lme4 1.1.37, BayesFactor 0.9.12-4.7) and Python 3.12 (numpy 1.26.4, pandas 2.2.2, scipy 1.17.1, matplotlib 3.10.7, statsmodels 0.14.6).

Run every command from the repository root. The first R script that sources the engine compiles `tddm_split.cpp` into `models/tDDM/.rcpp_cache/`, which takes about a minute once.

## 1. Data

The raw data are not redistributed here because the source repository carries no licence. Download and rebuild them:

```bash
python data/fetch_data.py
```

The script downloads `study4_ug.csv` and `study4_ec.csv` from [jpheffne/epe](https://github.com/jpheffne/epe) at a pinned commit. It then applies the original study's preregistered exclusion: a participant is kept when their rating of "neutral" falls inside the 100 × 100 px square at the centre of the emotion grid, which leaves 351 of 508. Each signal is divided by its SD over the retained trials. The SD is the population SD, and the signals are not centred. The script writes `data/study4_HC.csv` and `data/study4_dep.csv` and checks both against the SHA-256 hashes of the files the paper was fitted on. It stops if either hash differs.

The model reads `RPE_scale` and `EPE_X1_scale` (valence prediction error). Response times outside 0.15–10 s are dropped inside the engine (`tddm_load_study4`), which removes 281 of 7,020 trials.

## 2. Regenerate statistics and figures from the shipped fits (about 20 minutes)

```bash
bash run_fast.sh
```

This uses only the files in `results/` and runs no model fits. Logs go to `results/logs/`, figures to `docs/PAPER/figures/`. Section 4 lists what each step reproduces.

## 3. Refit from scratch (about 14 hours on 8 cores)

```bash
bash run_full.sh           # CORES=8 by default
```

Refits do not give bit-identical parameter estimates. The fitting scripts use `RNGkind("L'Ecuyer-CMRG")` with `mclapply`, so the random stream each participant receives depends on how participants are distributed over cores. DEoptim also resolves near-ties between candidate solutions differently when two likelihood implementations differ in the tenth decimal place. The group-level results are stable. Supplement S2 reports the spread across the ten runs of the main model, and Section 5 compares the shipped fits with a refit through the public fitting path.

## 4. Where each result comes from

`CANON` = `results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv`, the main model's per-participant estimates (mean signal first). For each participant it keeps the run with the lowest negative log-likelihood among 10 optimizer runs. The column `tau` holds the fitted value. Every analysis sets τ ≤ 0 to 0, because any τ ≤ 0 gives the same single-phase model.

Runtimes are wall-clock times on 8 to 12 cores (Apple M-series).

### Main model and baselines

| Paper item | Script and command | Output | Time |
|---|---|---|---|
| 10 optimizer runs of the main model, seeds 20260821–30 (the runs behind `CANON`) | `Rscript scripts/tDDM/diag_quadrature_refit.R --seed <seed> --cores 12`, once per seed | `results/tDDM/diag_quadrature_refit_20260927_*/fits_4j_accurate_rep01seed.csv` | ~26 min per run |
| `CANON` (best of 10) | `Rscript scripts/tDDM/memo_0919/build_accurate_bestrun.R` | `CANON`, `summary.txt` | s |
| 10 runs with the contrast signal first | `Rscript scripts/tDDM/optimizer_reproducibility.R --orientation contrast --reps 10 --cores 8` | `results/tDDM/optimizer_reproducibility_contrast_r1_20260928_005730/` | ~20 min |
| Contrast-first best of 10 | `Rscript scripts/tDDM/memo_0919/build_bestof_runs.R --runs-dir <that folder> --out <file>` | `.../canonical_4j_accurate_bestof10/fits_4j_accurate_contrast_first_best10.csv` | s |
| Baseline 4a (no contrast, no onset) | `Rscript scripts/tDDM/refit_corrected_likelihood.R --arms 4a_static1` | `results/tDDM/corrected_likelihood_20260819_012638/fits_4a_static1.csv` | <1 h |
| Contrast model 4c (no onset) | `Rscript scripts/tDDM/refit_corrected_likelihood.R --arms 4c_static` | `results/tDDM/corrected_likelihood_20260819_005659/fits_4c_static.csv` | <1 h |

The baselines have no onset, so their likelihood does not involve the split-at-τ integral described in Section 6.

### Results and Supplement

| Paper item | Script | Output |
|---|---|---|
| Group comparisons on drift, contrast weight, threshold, non-decision time, τ | `Rscript scripts/tDDM/group_stats.R --fits CANON` | stdout |
| Bound rates; τ ≤ 0 share and robustness; CES-D, age, gender; ANCOVA; BIC ladder (Fig. 3) | `Rscript scripts/tDDM/memo_0919/reanalysis_accurate.R` | `.../canonical_4j_accurate_bestof10/reanalysis_accurate.txt` |
| Model-free behaviour: acceptance, response time, trial-level choice model | `Rscript scripts/tDDM/memo_0919/model_free_bestrun.R` | stdout |
| Table 1 | `Rscript scripts/tDDM/demographics_table.R --fits CANON` | `results/tDDM/cohort_demographics/` |
| CES-D as a continuous score | `Rscript scripts/tDDM/ces_continuous_correlation.R` | stdout |
| S1 collinearity of the two signals | `Rscript scripts/paper/s1_collinearity.R` | `results/paper_numbers/` |
| Table S2 transition-density check | `Rscript scripts/tDDM/verify_transition_density.R` | `results/tDDM/transition_density_check_20260820_003321/` |
| S2 agreement among the 10 runs | `Rscript scripts/paper/s2_optimizer_runs_accurate.R` | stdout |
| Tables S3, S4, S6–S10 | `Rscript scripts/paper/supplement_tables_accurate.R` | `.../canonical_4j_accurate_bestof10/supplement_tables.txt` |
| Table S5, inversion of the observed Δτ, sensitivity ladder (Fig. 7) | `Rscript scripts/tDDM/memo_0919/stage6_summary_accurate.R` | `.../canonical_4j_accurate_bestof10/stage6_summary.txt` |
| Starting-point bound from acceptance rates (S4, Analysis 6) | `Rscript scripts/tDDM/memo_0919/starting_point_export.R && python scripts/tDDM/memo_0919/starting_point_sensitivity.py` | `results/tDDM/memo_analyses_0919/starting_point_sensitivity_20260927/` (~12 min) |

### Identifiability analyses (Supplement S4)

All were run by `scripts/tDDM/memo_0919/run_stage6_accurate_controls.sh` on 8 cores. Its log is `results/tDDM/stage6_accurate/runner.log`. Each simulated participant is fitted three times and the lowest-NLL fit is kept.

| Analysis | Script | Output | Time |
|---|---|---|---|
| 1. Individual recovery, 20 and 80 trials (Table S4, Fig. S1) | `recovery_redesign_tau.R --trans image --arms redesign --trials 20,80 --N 50 --reps 3` | `recovery_redesign_tau_20260928_192519/` | 0.4 h |
| 2. Split-half replication, mean-first and baseline (Fig. 5) | `split_half_predictive.R --K 30 --models 4a,4j_mean_first --trans image --itermax 100 --np 60 --restarts 2` | `split_half_predictive_image_20260928_115546/` | 1.9 h |
| 2. Split-half replication, contrast-first | the same with `--models 4j_contrast_first` | `split_half_predictive_image_20260928_134940/` | 0.7 h |
| 3. Positive control, 20 trials, 4 levels (Table S7, Fig. 7) | `positive_control_tau.R --reps 10 --gaps 0,0.085,0.169,0.34 --trans image --tau-base 0.804` | `positive_control_tau_20260928_231637/` | 2.2 h |
| 3. Positive control, 80 trials | the same with `--reps 5 --gaps 0,0.34 --ntrials 80` | `positive_control_tau_20260929_012719/` | 2.9 h |
| 4. Negative control, true Δτ = 0 (Tables S5–S6, Fig. 6) | `nulltau_corrected.R --reps 8 --tau-common 0.891` | `nulltau_corrected_image_20260928_211818/` | 1.1 h |
| 5. Response-time signature (Table S8, Fig. S2) | `rt_signature_model_check.R --fits CANON` | `rt_signature_check_20260928_005349/` | minutes |
| 6. Starting-point sensitivity, 0.05 and 0.10 | `nulltau_corrected.R --reps 5 --tau-common 0.891 --wgap 0.05` (and `0.10`) | `nulltau_corrected_image_20260928_222147/`, `..._224927/` | 0.5 h each |
| Non-decision-time bounds 1.5–3.5 s (Table S3) | `nDT_bound_ladder_engine.R --bounds 1.5,2.0,2.5,3.0,3.5` | `nDT_bound_ladder_engine_20260928_192521/` | 0.5 h |
| Choice signature | `memo_0919/choice_signature_check.R` | `choice_signature_check_20260928_005021/` | minutes |
| Effect of τ on the median response time | `tau_causes_rt_check.R` | `tau_causes_rt_20260928_005729/` | minutes |

`--tau-common 0.891` and `--tau-base 0.804` are the pooled and the healthy-control mean onsets in `CANON` with τ ≤ 0 set to 0.

### Figures

All figure scripts write to `docs/PAPER/figures/`.

| Figure | Script |
|---|---|
| 1, 2 (schematics, no data) | `python scripts/visualization/make_schematic_figures.py` |
| 3, 4, 5 | `python scripts/visualization/make_paper_figures_v2.py` |
| 6, 7 | `python scripts/visualization/make_boundary_figures.py` |
| S1 | `python scripts/tDDM/figS1_individual_recovery.py` |
| S2 | `python scripts/tDDM/figS2_absolute_fit.py` (run `rt_signature_model_check.R` first; its trial-level file contains observed response times and is not shipped) |
| S3, S4 | `python scripts/eeg/figs_cesd_and_age.py` |

## 5. Reproducing the main fits through the public fitting path

`CANON` was fitted with `scripts/tDDM/diag_quadrature_refit.R`, which evaluates the time-varying likelihood in R. `run_full.sh` fits through the engine instead (`optimizer_reproducibility.R`), which evaluates the same likelihood in C++ and is 10 to 17 times faster. The two implementations agree to within 1.3 × 10⁻¹⁰ in log-likelihood on the stored fits, on 3,000 random parameter points, and on simulated data (`verify_cpp_likelihood.R`, output in `canonical_4j_accurate_bestof10/verify_cpp_likelihood.txt`).

We refitted the main model 10 times through the engine (`results/tDDM/optimizer_reproducibility_mean_r1_20261002_130237/`, 17 minutes on 12 cores) and took the best run per participant with `build_bestof_runs.R`:

| | Δτ (DEP − HC) | p | mean ΔBIC (4j − 4c) | ΔBIC < 0 | ΔBIC < −2 |
|---|---|---|---|---|---|
| `CANON` (paper) | +0.209 s | .044 | +1.790 | 10.0% | 2.8% |
| engine refit | +0.226 s | .031 | +1.789 | 10.0% | 2.8% |

The best negative log-likelihood agreed to within 10⁻⁶ for 84.9% of participants and to within 0.01 for all but two, where the engine run found a lower value (by up to 0.125). Fitted onsets correlated at r = 0.983. The engine refit's Δτ lies inside the range the paper reports across its own ten runs (+0.171 to +0.280 s, Supplement S2).

## 6. Notes on the likelihood

- `models/tDDM/tddm_engine.R` is the only implementation of the likelihood. The time-varying likelihood integrates the state density at τ (method of images) times the second-phase first-passage density over the evidence level. The engine does this by adaptive Gauss–Kronrod quadrature (relative and absolute tolerance 10⁻⁶) with breakpoints at a distance of min(bsep/4, 20√(t − τ)) from each boundary, because the second-phase density concentrates near the boundaries when a response falls shortly after τ. Supplement S1 describes the method, and `diag_quadrature_spectral.R`, `diag_quadrature_brute.R`, `diag_quadrature_validate.R`, and `diag_quadrature_mc_truth.py` check it against a spectral solution, a brute-force integral, and Monte Carlo simulation.
- `TDDM_QUAD` in the engine also accepts `"grid"`, a fixed 30-point rule we used earlier. It overstates the likelihood of a response that falls just after τ, and no reported result uses it. It is kept so that the comparison in `diag_quadrature_refit_*/summary.txt` can be reproduced.
- The engine's `--trans pipeline` option keeps an earlier transition density that ignored absorption and mis-scaled the diffusion. It exists so that the density comparison in Table S2 can be reproduced. Every reported fit uses `--trans image`.
- The trial simulator used by the identifiability analyses (`tddm_simulate_trial`, `tddm_simulate_block`) steps in 1 ms, ends a path with the Brownian-bridge probability that it crossed a boundary within the step, and cuts the step at τ. `diag_simulator_vs_likelihood.R` checks its decision-time and choice distributions against the likelihood.
- DEoptim runs with its default value-to-reach, so every fit runs all 200 generations. Setting `VTR` in `DEoptim.control()` makes DEoptim stop at generation 0.

## Citation

Please cite the paper and the original dataset:

Heffner, J., Son, J.-Y., & FeldmanHall, O. (2021). Emotion prediction errors guide socially adaptive behaviour. *Nature Human Behaviour*, 5, 1391–1401.

## Licence

Code: MIT (see `LICENSE`). The derived per-participant estimates and simulation outputs in `results/` are released under CC BY 4.0. The behavioural data belong to the original authors and are not included.
