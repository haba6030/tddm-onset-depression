# config.R — Source this from any script
# Usage: source("../../config.R") from scripts/tDDM/, or source("config.R") from root

PROJECT_ROOT <- tryCatch(
  normalizePath(file.path(dirname(sys.frame(1)$ofile), "."), mustWork = FALSE),
  error = function(e) ""
)

# Walk up from script location to find config.R at project root
if (!nzchar(PROJECT_ROOT) || PROJECT_ROOT == ".") {
  # Try here package
  if (requireNamespace("here", quietly = TRUE)) {
    PROJECT_ROOT <- here::here()
  } else {
    # Fallback: search parent directories for config.R
    candidate <- getwd()
    while (nchar(candidate) > 1) {
      if (file.exists(file.path(candidate, "config.R")) &&
          file.exists(file.path(candidate, "HKUSRP.Rproj"))) {
        PROJECT_ROOT <- candidate
        break
      }
      candidate <- dirname(candidate)
    }
  }
}

DATA_DIR    <- file.path(PROJECT_ROOT, "data")
MODEL_DIR   <- file.path(PROJECT_ROOT, "models")
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
FIGURES_DIR <- file.path(PROJECT_ROOT, "figures")
SCRIPTS_DIR <- file.path(PROJECT_ROOT, "scripts")
