# =============================================================================
# 04_robustness_spj.R  —  SPJ across four subsamples  (Table 3, robustness)
# -----------------------------------------------------------------------------
# Re-estimates the SPJ-corrected RTA elasticity on four subsamples to confirm
# the sign and rough magnitude are stable:
#   Pre-2008 (1995-2007), Post-2008 (2008-2020), Non-EU pairs, WTO members only.
# These form the robustness rows of Table 3.
# =============================================================================

source("code/01_packages_data.R")
dir.create("output", showWarnings = FALSE)

# SPJ on an arbitrary subsample, returning coef/se/stars for all regressors.
run_spj_sub <- function(dat, regs, label, seed = 42) {
  set.seed(seed)
  cat(sprintf("  %-24s (%s obs)...", label, format(nrow(dat), big.mark = ",")))
  ctrs <- unique(dat$exp_iso3)
  if (length(ctrs) < 6) { cat(" SKIP (too few countries)\n"); return(NULL) }

  fit_f <- tryCatch(
    fepois(as.formula(paste("trade_kusd ~", paste(regs, collapse = "+"),
                            "| pair_id + exp_time_id + imp_time_id")),
           data = dat, vcov = ~pair_id, glm.iter = 200, nthreads = 1),
    error = function(e) NULL)
  if (is.null(fit_f)) { cat(" FAIL\n"); return(NULL) }

  b_f  <- coef(fit_f)[regs]
  ctrs <- sample(ctrs); N <- length(ctrs)
  gA   <- ctrs[seq_len(floor(N / 2))]
  gB   <- ctrs[seq_len(N) > floor(N / 2)]
  b_ss <- lapply(list(AA = list(e = gA, i = gA), AB = list(e = gA, i = gB),
                       BA = list(e = gB, i = gA), BB = list(e = gB, i = gB)),
                 function(sp) {
    sd <- dat %>%
      filter(exp_iso3 %in% sp$e, imp_iso3 %in% sp$i) %>%
      mutate(p_id = as.integer(factor(paste(exp_iso3, imp_iso3))),
             e_id = as.integer(factor(paste(exp_iso3, year))),
             i_id = as.integer(factor(paste(imp_iso3, year))))
    fs <- tryCatch(
      fepois(as.formula(paste("trade_kusd ~", paste(regs, collapse = "+"),
                              "| p_id + e_id + i_id")),
             data = sd, vcov = ~p_id, glm.iter = 200, nthreads = 1),
      error = function(e) NULL)
    if (!is.null(fs)) coef(fs)[regs] else NULL
  })
  b_ss <- Filter(Negate(is.null), b_ss)
  if (length(b_ss) < 2) { cat(" FAIL (too few subpanels)\n"); return(NULL) }

  b_spj      <- 2 * b_f - Reduce("+", b_ss) / length(b_ss)
  available_r <- intersect(regs, names(b_f))
  vc <- tryCatch(
    vcov(fit_f, type = "cluster", cluster = ~pair_id)[available_r, available_r, drop = FALSE],
    error = function(e) diag(length(available_r)) * NA_real_)

  se <- rep(NA_real_, length(regs)); names(se) <- regs
  se[available_r] <- sqrt(diag(vc))
  pv <- rep(NA_real_, length(regs)); names(pv) <- regs
  pv[available_r] <- 2 * pnorm(-abs(b_spj[available_r] / se[available_r]))

  cat(" OK\n")
  tibble(sample = label, n_obs = nrow(dat), variable = names(b_spj),
         coef = round(b_spj, 6), se = round(se, 6), pval = round(pv, 6),
         stars = add_stars(pv))
}

rob_list <- list()
rob_list[["pre2008"]]  <- run_spj_sub(panel %>% filter(year <= 2007),
                                      REGRESSORS, "Pre-2008 (1995-2007)")
rob_list[["post2008"]] <- run_spj_sub(panel %>% filter(year >= 2008),
                                      REGRESSORS, "Post-2008 (2008-2020)")
rob_list[["non_eu"]]   <- run_spj_sub(panel %>% filter(!(eu_o == 1 & eu_d == 1)),
                                      REGRESSORS, "Non-EU pairs")
rob_list[["wto_only"]] <- run_spj_sub(panel %>% filter(wto_o == 1, wto_d == 1),
                                      REGRESSORS, "WTO members only")

results_rob <- bind_rows(rob_list)

if (nrow(results_rob) > 0) {
  t3_rob <- results_rob %>%
    mutate(cell = mapply(fmt_cell, coef, se, stars),
           Variable = dplyr::case_when(
             variable == "rta"         ~ "RTA",
             variable == "ln_dist"     ~ "ln(distance)",
             variable == "contig"      ~ "Contiguity",
             variable == "comlang_off" ~ "Common language",
             variable == "comcol"      ~ "Common coloniser",
             TRUE ~ variable)) %>%
    select(Variable, sample, cell) %>%
    distinct(Variable, sample, .keep_all = TRUE) %>%
    pivot_wider(names_from = sample, values_from = cell, values_fn = first)

  write.csv(results_rob, "output/table3_robustness_long.csv", row.names = FALSE)
  write.csv(t3_rob,      "output/table3_robustness_wide.csv", row.names = FALSE)
}

cat("\n04_robustness_spj.R complete: wrote Table 3 robustness block.\n\n")
