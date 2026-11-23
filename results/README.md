# results

`results/tDDM/` holds the per-participant estimates and the output folders of the runs the paper reports. The top-level README (Section 4) lists which folder supports which paper item and the command that wrote it.

A run folder is named `<script>_<YYYYMMDD_HHMMSS>`, after the time the run started. The three runs added in November 2026 (`recovery_baseline_20261106_210914`, `split_half_predictive_image_20261106_205634`, `tau_causes_rt_clampfirst_20261107_104449`) also contain `run_info.txt` with the command, commit and package versions. For the earlier runs the command is given in the top-level README.

| Folder | Contents |
|---|---|
| `canonical_4j_accurate_bestof10/` | `CANON` and the contrast-first best of 10, plus the summary files the paper quotes: `summary.txt`, `baseline_primary.txt`, `reanalysis_accurate.txt`, `stage6_summary.txt`, `supplement_tables.txt`, `verify_cpp_likelihood.txt` |
| `corrected_likelihood_20260819_*/` | baseline (`fits_4a_static1.csv`) and contrast (`fits_4c_static.csv`) model fits |
| `diag_quadrature_refit_20260927_*/` | the ten runs of the time-varying model behind `CANON` (seeds 20260821–30) |
| `optimizer_reproducibility_*/` | ten engine runs per signal ordering (contrast first; mean first for the check in Section 5) |
| `recovery_*`, `split_half_*`, `nulltau_*`, `positive_control_*`, `rt_signature_*`, `choice_signature_*`, `tau_causes_rt_*`, `nDT_bound_ladder_*` | identifiability and validation analyses (Supplement S3, S4) |
| `transition_density_check_*`, `diag_quadrature_spectral_*` | checks of the likelihood (Table S2; Section 6 of the top-level README) |

Two folders the paper does not use were removed on 2026-11-23 and remain in commit `ed36e00`: `canonical_4j_accurate_bestof10/clusters/` (a clustering of the fitted parameters) and `tau_causes_rt_20260928_005729/` (the onset-shift check before onsets at or below zero were set to zero; the paper reports `tau_causes_rt_clampfirst_20261107_104449`).

Files that contain observed trial-level response times (`trial_level.csv`, `trials.csv`) are not shipped. Running the scripts writes further folders (`results/logs/`, `results/paper_numbers/`, `results/tDDM/cohort_demographics/`, new timestamped run folders), which are not tracked.
