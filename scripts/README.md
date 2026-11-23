# scripts

Every script runs from the repository root and reads its paths from `config.R`. The top-level README (Section 4) maps each paper item to a script and gives the command.

The folder and file names come from the development repository the paper was written in. They are kept unchanged so that the paths match the run records and the summary files in `results/`.

| Folder | What the scripts do |
|---|---|
| `tDDM/` | Model fits: `refit_corrected_likelihood.R` (baseline and contrast models), `diag_quadrature_refit.R` and `optimizer_reproducibility.R` (time-varying model). Group statistics: `group_stats.R`, `demographics_table.R`, `ces_continuous_correlation.R`. Identifiability and validation analyses: `recovery_baseline.R`, `recovery_redesign_tau.R`, `split_half_predictive.R`, `nulltau_corrected.R`, `positive_control_tau.R`, `rt_signature_model_check.R`, `tau_causes_rt_check.R`, `nDT_bound_ladder_engine.R`. Checks of the likelihood and the simulator: `verify_*.R`, `diag_quadrature_spectral.R`, `diag_simulator_vs_likelihood.R`. Figures S1 and S4: `figS1_individual_recovery.py`, `figS2_absolute_fit.py` |
| `tDDM/memo_0919/` | Best-of-10 selection (`build_accurate_bestrun.R`, `build_bestof_runs.R`), reanalysis and summary of the time-varying fits (`reanalysis_accurate.R`, `stage6_summary_accurate.R`), model-free comparison, choice signature, starting-point calibration. The folder is named after the analysis memo of 19 September 2026 that introduced these scripts |
| `paper/` | Numbers for the baseline model (`baseline_primary_accurate.R`) and for the Supplement (`s1_collinearity.R`, `s2_optimizer_runs_accurate.R`, `supplement_tables_accurate.R`) |
| `visualization/` | Figs. 1–4, S2 and S3 |
| `eeg/` | Figs. S5 and S6 (CES-D as a continuous measure, age). The folder name is historical; these figures use behavioural fits only |

The likelihood itself lives in `models/tDDM/tddm_engine.R`, and every fitting script sources it from there.
