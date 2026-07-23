# =============================================================================
# 02_descriptives.R  —  Summary statistics and fixed-effect counts
# -----------------------------------------------------------------------------
# Reproduces the descriptive figures cited in Section 4.1 (panel size, zero
# share, variance-to-mean ratio, fixed-effect level counts). Output feeds the
# in-text data description; not a numbered results table.
# =============================================================================

source("code/01_packages_data.R")
dir.create("output", showWarnings = FALSE)

# ---- Summary statistics ----------------------------------------------------
sumstats <- tibble(
  Variable = c("Trade (USD thous.) — all", "Trade (USD thous.) — non-zero",
               "RTA", "ln(distance)", "Contiguity", "Common language",
               "Common coloniser", "EU exporter", "WTO exporter",
               "Var/mean ratio (non-zero)"),
  N    = c(nrow(panel), length(nz), rep(nrow(panel), 7), length(nz)),
  Mean = round(c(mean(panel$trade_kusd), mean(nz),
                 mean(panel$rta), mean(panel$ln_dist), mean(panel$contig),
                 mean(panel$comlang_off), mean(panel$comcol),
                 mean(panel$eu_o), mean(panel$wto_o), var(nz) / mean(nz)), 3),
  SD   = round(c(sd(panel$trade_kusd), sd(nz),
                 sd(panel$rta), sd(panel$ln_dist), sd(panel$contig),
                 sd(panel$comlang_off), sd(panel$comcol),
                 sd(panel$eu_o), sd(panel$wto_o), NA_real_), 3),
  Min  = round(c(min(panel$trade_kusd), min(nz), 0, min(panel$ln_dist),
                 0, 0, 0, 0, 0, NA_real_), 3),
  Max  = round(c(max(panel$trade_kusd), max(nz), 1, max(panel$ln_dist),
                 1, 1, 1, 1, 1, NA_real_), 3)
)

# ---- Fixed-effect level counts --------------------------------------------
fe_tbl <- tibble(
  `Fixed effect`  = c("Pair (mu_ij)", "Exporter-time (lambda_it)",
                      "Importer-time (xi_jt)", "Year"),
  `Column`        = c("pair_id", "exp_time_id", "imp_time_id", "year"),
  `Unique levels` = c(n_distinct(panel$pair_id), n_distinct(panel$exp_time_id),
                      n_distinct(panel$imp_time_id), n_distinct(panel$year))
)

write.csv(sumstats, "output/table_descriptives_summary.csv", row.names = FALSE)
write.csv(fe_tbl,   "output/table_descriptives_fe_counts.csv", row.names = FALSE)

cat("02_descriptives.R complete: wrote summary statistics and FE counts.\n\n")
