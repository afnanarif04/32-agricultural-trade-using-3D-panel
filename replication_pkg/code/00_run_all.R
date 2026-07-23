# =============================================================================
# 00_run_all.R  —  Master replication script
# -----------------------------------------------------------------------------
# Reproduces every table in the paper. Set the working directory to the
# REPOSITORY ROOT (the folder containing code/, data/, output/) before running,
# then source this file. Each script is self-contained and also sources
# 01_packages_data.R on its own, so scripts can be run individually too.
#
#   Option 1 (this file): source("code/00_run_all.R")
#   Option 2 (per table) : source("code/03_main_estimation.R"), etc.
# =============================================================================

rm(list = ls())
set.seed(42)

code_dir <- "code"

scripts <- c(
  "01_packages_data.R",     # load panel, build objects, shared helpers
  "02_descriptives.R",      # summary statistics, FE counts  (paper Section 4.1)
  "03_main_estimation.R",   # nine estimators on full panel  (Table 3, full panel)
  "04_robustness_spj.R",    # SPJ across four subsamples      (Table 3, robustness)
  "05_monte_carlo.R"        # simulation over 81 DGP cells    (Table 1)
)

cat("=====================================================================\n")
cat("  Replication master runner\n")
cat("=====================================================================\n\n")

for (s in scripts) {
  cat(sprintf("--- Running %-28s", s))
  t0 <- proc.time()["elapsed"]
  tryCatch(
    source(file.path(code_dir, s), echo = FALSE),
    error = function(e) cat("\n  *** ERROR:", conditionMessage(e), "\n")
  )
  cat(sprintf(" [%.1f s]\n", proc.time()["elapsed"] - t0))
}

cat("\n=====================================================================\n")
cat("  All scripts complete. Output tables written to output/.\n")
cat("=====================================================================\n")
