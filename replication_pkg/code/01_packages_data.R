# =============================================================================
# 01_packages_data.R  —  Packages, data loading, shared objects and helpers
# -----------------------------------------------------------------------------
# Sourced first by every downstream script. Loads the merged bilateral
# agricultural trade panel, applies the analysis-ready renames, and defines the
# formulas, regressor list, and the star/format helpers used everywhere.
# =============================================================================

# ---- 1. Packages: auto-install if missing ---------------------------------
# (All result tables are written as CSV; openxlsx is not required.)
required_pkgs <- c("data.table", "dplyr", "tidyr", "fixest")
to_install <- required_pkgs[!vapply(required_pkgs,
                                    requireNamespace, logical(1),
                                    quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(fixest)
})

# ---- 2. Masking guards (re-bind the dplyr verbs) --------------------------
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise

# ---- 3. Single-thread guard (hardware-safe, reproducible) -----------------
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1)
}

# ---- 4. Load and validate the panel ---------------------------------------
panel <- fread("data/BiasCorr3DPPML_Agri_MergedPanel.csv") %>%
  as_tibble() %>%
  rename(
    exp_iso3   = iso3_o,
    imp_iso3   = iso3_d,
    trade_kusd = y_ijt,
    ln_trade   = ln_y_ijt,
    rta        = fta_wto
  )

stopifnot(nrow(panel) == 40560)
stopifnot(all(c("exp_iso3", "imp_iso3", "year", "trade_kusd", "ln_trade",
                "pair_id", "exp_time_id", "imp_time_id", "rta",
                "ln_dist", "contig", "comlang_off", "comcol") %in% names(panel)))

# non-zero trade vector, reused for dispersion diagnostics
nz <- panel$trade_kusd[panel$trade_kusd > 0]

# ---- 5. Shared model objects ----------------------------------------------
# The paper's five reported gravity regressors. Time-invariant dyadic controls
# (contig, comlang_off, comcol) are absorbed by the pair fixed effect in the
# three-way specifications and identified only in the two-way comparator.
REGRESSORS <- c("rta", "ln_dist", "contig", "comlang_off", "comcol")

fml_3d  <- as.formula(paste("trade_kusd ~", paste(REGRESSORS, collapse = "+"),
                            "| pair_id + exp_time_id + imp_time_id"))
fml_2w  <- as.formula(paste("trade_kusd ~", paste(REGRESSORS, collapse = "+"),
                            "| pair_id + year"))
fml_ols <- as.formula(paste("ln_trade ~",   paste(REGRESSORS, collapse = "+"),
                            "| pair_id + exp_time_id + imp_time_id"))

# Reporting order of the nine estimators
EST_ORDER <- c("SPJ", "3D-PPML", "WZ-Analytical", "Polyad", "Bootstrap",
               "Heckman-2S", "GNLS", "OLS-3D", "2W-PPML")

# ---- 6. Shared helpers -----------------------------------------------------
# Significance stars from a p-value.
add_stars <- function(pv) {
  dplyr::case_when(pv < 0.01 ~ "***", pv < 0.05 ~ "**",
                   pv < 0.10 ~ "*",  TRUE ~ "")
}

# Extract coef / se / p-value tibble from a fixest fit, robust to dropped
# (collinear or singleton-absorbed) regressors.
get_results <- function(fit, label, regs = REGRESSORS) {
  if (is.null(fit)) return(NULL)
  tryCatch({
    available <- intersect(regs, names(coef(fit)))
    if (length(available) == 0) {
      message(sprintf("[%s] all regressors dropped", label))
      return(NULL)
    }
    cf <- coef(fit)[available]
    vc <- tryCatch(
      vcov(fit, type = "cluster", cluster = ~pair_id)[available, available, drop = FALSE],
      error = function(e) vcov(fit)[available, available, drop = FALSE]
    )
    se <- sqrt(diag(vc))
    pv <- 2 * pnorm(-abs(cf / se))
    out <- tibble(estimator = label, variable = regs,
                  coef = NA_real_, se = NA_real_, pval = NA_real_, stars = "")
    idx <- match(available, regs)
    out$coef[idx]  <- round(cf, 6)
    out$se[idx]    <- round(se, 6)
    out$pval[idx]  <- round(pv, 6)
    out$stars[idx] <- add_stars(pv)
    out
  }, error = function(e) {
    message(sprintf("[%s] %s", label, e$message)); NULL
  })
}

# "coef*** (se)" cell formatter for the wide result tables.
fmt_cell <- function(coef, se, stars) {
  if (is.na(coef)) return("--")
  sprintf("%.3f%s (%.3f)", coef, stars, se)
}

cat("01_packages_data.R complete: panel loaded (",
    format(nrow(panel), big.mark = ","), " obs).\n\n", sep = "")
