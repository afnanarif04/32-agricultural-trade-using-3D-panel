# =============================================================================
# 03_main_estimation.R  —  Nine estimators on the full panel  (Table 3, full)
# -----------------------------------------------------------------------------
# Estimates the regional-trade-agreement (RTA) elasticity on the full panel of
# 40,560 observations with nine estimators:
#   1 SPJ            split-panel-jackknife bias correction (proposed)
#   2 3D-PPML        uncorrected three-way Poisson PMLE
#   3 WZ-Analytical  Weidner-Zylkin (2021) leave-one-pair-out approximation
#   4 Polyad         Resende et al. (2025) — placeholder (code not public)
#   5 Bootstrap      Higgins & Jochmans (2024) pair bootstrap correction
#   6 Heckman-2S     Mnasri & Nechi (2021) selection two-step
#   7 GNLS           Mnasri & Nechi (2021) iterated GNLS
#   8 OLS-3D         log-linearised OLS with three-way FE
#   9 2W-PPML        two-way Poisson (pair + year) comparator
#
# Five of these (SPJ, 3D-PPML, WZ-Analytical, GNLS, OLS-3D) are reported as
# columns in Table 3; the rest are discussed in Section 4.5.
#
# NOTE ON RUNTIME: full bootstrap uses B = 200 and the run takes noticeably
# longer than the other estimators. Set B_BOOT below to a smaller value for a
# quick check. All results in the paper use B_BOOT = 200.
# =============================================================================

source("code/01_packages_data.R")
dir.create("output", showWarnings = FALSE)

set.seed(42)
B_BOOT <- 200

results <- list()

# ---- [2] 3D-PPML (fitted first; SPJ, WZ, Bootstrap, GNLS depend on it) -----
cat("  [1/9] 3D-PPML...\n")
fit_ppml <- tryCatch(
  fepois(fml_3d, data = panel, vcov = ~pair_id, glm.iter = 200, nthreads = 1),
  error = function(e) { cat("    FAIL:", e$message, "\n"); NULL }
)
results[["PPML_3D"]] <- get_results(fit_ppml, "3D-PPML")

# ---- [1] SPJ: split-panel jackknife ---------------------------------------
# Split the country set into two halves along the cross-sectional dimension;
# the four exporter-half x importer-half subpanels each carry twice the
# leading bias, so b_spj = 2*b_full - mean(b_subpanels) cancels it to O(1/N^2).
cat("  [2/9] SPJ...\n")
run_spj <- function(dat, fml_full, regs, seed = 42) {
  set.seed(seed)
  fit_full <- tryCatch(
    fepois(fml_full, data = dat, vcov = ~pair_id, glm.iter = 200, nthreads = 1),
    error = function(e) NULL
  )
  if (is.null(fit_full)) return(NULL)
  b_full <- coef(fit_full)[regs]
  ctrs   <- sample(unique(dat$exp_iso3))
  N      <- length(ctrs)
  gA     <- ctrs[seq_len(floor(N / 2))]
  gB     <- ctrs[seq_len(N) > floor(N / 2)]
  subs   <- list(AA = list(e = gA, i = gA), AB = list(e = gA, i = gB),
                 BA = list(e = gB, i = gA), BB = list(e = gB, i = gB))
  b_subs <- lapply(subs, function(sp) {
    sd <- dat %>%
      filter(exp_iso3 %in% sp$e, imp_iso3 %in% sp$i) %>%
      mutate(p_id = as.integer(factor(paste(exp_iso3, imp_iso3))),
             e_id = as.integer(factor(paste(exp_iso3, year))),
             i_id = as.integer(factor(paste(imp_iso3, year))))
    fm <- as.formula(paste("trade_kusd ~", paste(regs, collapse = "+"),
                           "| p_id + e_id + i_id"))
    fs <- tryCatch(
      fepois(fm, data = sd, vcov = ~p_id, glm.iter = 200, nthreads = 1),
      error = function(e) NULL)
    if (!is.null(fs)) coef(fs)[regs] else NULL
  })
  b_subs <- Filter(Negate(is.null), b_subs)
  if (length(b_subs) < 2) return(NULL)
  b_avg <- Reduce("+", b_subs) / length(b_subs)
  b_spj <- 2 * b_full - b_avg
  vc <- tryCatch(
    vcov(fit_full, type = "cluster", cluster = ~pair_id)[regs, regs, drop = FALSE],
    error = function(e) vcov(fit_full)[regs, regs, drop = FALSE])
  list(coef = b_spj, vcov = vc, n_sub = length(b_subs))
}
spj <- run_spj(panel, fml_3d, REGRESSORS)
if (!is.null(spj)) {
  cf <- spj$coef; se <- sqrt(diag(spj$vcov)); pv <- 2 * pnorm(-abs(cf / se))
  results[["SPJ"]] <- tibble(estimator = "SPJ", variable = names(cf),
                             coef = round(cf, 6), se = round(se, 6),
                             pval = round(pv, 6), stars = add_stars(pv))
}

# ---- [3] WZ-Analytical: leave-one-pair-out approximation ------------------
# The score-based WZ correction is degenerate at the PPML optimum (the score is
# numerically zero), so a leave-one-pair-out jackknife is used as a stable,
# consistent approximation to the WZ (2021) correction.
cat("  [3/9] WZ-Analytical...\n")
run_wz <- function(fit_full, dat, regs) {
  if (is.null(fit_full)) return(NULL)
  tryCatch({
    b_full <- coef(fit_full)[regs]
    avail  <- intersect(regs, names(b_full))
    set.seed(99)
    pairs_sample <- sample(unique(dat$pair_id),
                           min(40, n_distinct(dat$pair_id)))
    bias_list <- lapply(pairs_sample, function(pid) {
      dat_sub <- dat %>% filter(pair_id != pid)
      fml_s <- as.formula(paste("trade_kusd ~", paste(regs, collapse = "+"),
                                "| pair_id + exp_time_id + imp_time_id"))
      fit_s <- tryCatch(
        fepois(fml_s, data = dat_sub, vcov = ~pair_id, glm.iter = 100, nthreads = 1),
        error = function(e) NULL)
      if (is.null(fit_s)) return(NULL)
      avail_s <- intersect(avail, names(coef(fit_s)))
      b_s <- rep(NA_real_, length(avail)); names(b_s) <- avail
      b_s[avail_s] <- coef(fit_s)[avail_s]
      b_s - b_full[avail]
    })
    valid_lopo <- Filter(Negate(is.null), bias_list)
    vc <- tryCatch(
      vcov(fit_full, type = "cluster", cluster = ~pair_id)[regs, regs, drop = FALSE],
      error = function(e) vcov(fit_full)[regs, regs, drop = FALSE])
    if (length(valid_lopo) < 5) return(list(coef = b_full, vcov = vc))
    N_pairs   <- n_distinct(dat$pair_id)
    bias_mean <- Reduce("+", valid_lopo) / length(valid_lopo) * (N_pairs - 1)
    b_wz      <- b_full
    b_wz[avail] <- b_full[avail] - bias_mean
    list(coef = b_wz, vcov = vc)
  }, error = function(e) { cat("    WZ error:", e$message, "\n"); NULL })
}
wz <- run_wz(fit_ppml, panel, REGRESSORS)
if (!is.null(wz)) {
  cf <- wz$coef; se <- sqrt(diag(wz$vcov)); pv <- 2 * pnorm(-abs(cf / se))
  results[["WZ"]] <- tibble(estimator = "WZ-Analytical", variable = names(cf),
                            coef = round(cf, 6), se = round(se, 6),
                            pval = round(pv, 6), stars = add_stars(pv))
}

# ---- [4] Polyad placeholder (Resende et al. 2025; code not public) --------
cat("  [4/9] Polyad: placeholder (public code unavailable)\n")
results[["Polyad"]] <- tibble(estimator = "Polyad", variable = REGRESSORS,
                              coef = NA_real_, se = NA_real_,
                              pval = NA_real_, stars = "")

# ---- [5] Bootstrap bias correction (Higgins & Jochmans 2024) --------------
cat(sprintf("  [5/9] Bootstrap (B = %d)...\n", B_BOOT))
run_boot <- function(fit_full, dat, regs, B = 200, seed = 42) {
  if (is.null(fit_full)) return(NULL)
  set.seed(seed)
  b_full <- coef(fit_full)[regs]
  pairs  <- unique(dat$pair_id)
  b_draws <- lapply(seq_len(B), function(b) {
    drawn <- sample(pairs, length(pairs), replace = TRUE)
    bdat  <- dat %>%
      filter(pair_id %in% drawn) %>%
      mutate(p_id = as.integer(factor(paste(exp_iso3, imp_iso3))),
             e_id = as.integer(factor(paste(exp_iso3, year))),
             i_id = as.integer(factor(paste(imp_iso3, year))))
    fm <- as.formula(paste("trade_kusd ~", paste(regs, collapse = "+"),
                           "| p_id + e_id + i_id"))
    fs <- tryCatch(
      fepois(fm, data = bdat, vcov = ~p_id, glm.iter = 100, nthreads = 1),
      error = function(e) NULL)
    if (!is.null(fs)) coef(fs)[regs] else NULL
  })
  valid <- Filter(Negate(is.null), b_draws)
  if (length(valid) < 2) return(NULL)
  b_bar  <- Reduce("+", valid) / length(valid)
  b_boot <- 2 * b_full - b_bar
  vc <- tryCatch(
    vcov(fit_full, type = "cluster", cluster = ~pair_id)[regs, regs, drop = FALSE],
    error = function(e) vcov(fit_full)[regs, regs, drop = FALSE])
  list(coef = b_boot, vcov = vc, n_draws = length(valid))
}
boot <- run_boot(fit_ppml, panel, REGRESSORS, B = B_BOOT)
if (!is.null(boot)) {
  cf <- boot$coef; se <- sqrt(diag(boot$vcov)); pv <- 2 * pnorm(-abs(cf / se))
  results[["Boot"]] <- tibble(estimator = "Bootstrap", variable = names(cf),
                              coef = round(cf, 6), se = round(se, 6),
                              pval = round(pv, 6), stars = add_stars(pv))
}

# ---- [6] Heckman two-step (Mnasri & Nechi 2021) ---------------------------
cat("  [6/9] Heckman-2S...\n")
run_heck <- function(dat, regs) {
  tryCatch({
    dat2 <- dat %>% mutate(trade_pos = as.integer(trade_kusd > 0))
    s1 <- glm(as.formula(paste("trade_pos ~", paste(regs, collapse = "+"))),
              data = dat2, family = binomial("probit"))
    xb <- predict(s1)
    dat2$imr <- dnorm(xb) / pmax(pnorm(xb), 1e-10)
    dat2$imr[!is.finite(dat2$imr)] <- 0
    fml_s2 <- as.formula(paste("trade_kusd ~", paste(c(regs, "imr"), collapse = "+"),
                               "| pair_id + exp_time_id + imp_time_id"))
    s2 <- fepois(fml_s2, data = dat2, vcov = ~pair_id, glm.iter = 200, nthreads = 1)
    vc <- tryCatch(
      vcov(s2, type = "cluster", cluster = ~pair_id)[regs, regs, drop = FALSE],
      error = function(e) vcov(s2)[regs, regs, drop = FALSE])
    list(coef = coef(s2)[regs], vcov = vc)
  }, error = function(e) { cat("    Heckman error:", e$message, "\n"); NULL })
}
heck <- run_heck(panel, REGRESSORS)
if (!is.null(heck)) {
  cf <- heck$coef; se <- sqrt(diag(heck$vcov)); pv <- 2 * pnorm(-abs(cf / se))
  results[["Heckman2S"]] <- tibble(estimator = "Heckman-2S", variable = names(cf),
                                   coef = round(cf, 6), se = round(se, 6),
                                   pval = round(pv, 6), stars = add_stars(pv))
}

# ---- [7] Iterated GNLS (Mnasri & Nechi 2021) ------------------------------
cat("  [7/9] GNLS...\n")
run_gnls <- function(fit_init, dat, fml_full, regs, max_iter = 10, tol = 1e-6) {
  if (is.null(fit_init)) return(NULL)
  tryCatch({
    mu_full <- predict(fit_init, newdata = dat)
    mu_full[is.na(mu_full)] <- mean(fitted(fit_init), na.rm = TRUE)
    b_p   <- coef(fit_init)[intersect(regs, names(coef(fit_init)))]
    fit_c <- fit_init
    for (k in seq_len(max_iter)) {
      dat_w <- dat %>% mutate(gnls_wt = 1 / pmax(mu_full, 1e-8))
      fit_k <- tryCatch(
        fepois(fml_full, data = dat_w, weights = ~gnls_wt,
               vcov = ~pair_id, glm.iter = 200, nthreads = 1),
        error = function(e) NULL)
      if (is.null(fit_k)) break
      b_k <- coef(fit_k)[intersect(regs, names(coef(fit_k)))]
      if (max(abs(b_k - b_p[names(b_k)]), na.rm = TRUE) < tol) { fit_c <- fit_k; break }
      mu_full <- predict(fit_k, newdata = dat)
      mu_full[is.na(mu_full)] <- mean(fitted(fit_k), na.rm = TRUE)
      b_p <- b_k; fit_c <- fit_k
    }
    avail <- intersect(regs, names(coef(fit_c)))
    vc <- tryCatch(
      vcov(fit_c, type = "cluster", cluster = ~pair_id)[avail, avail, drop = FALSE],
      error = function(e) vcov(fit_c)[avail, avail, drop = FALSE])
    cf_out <- rep(NA_real_, length(regs)); names(cf_out) <- regs
    cf_out[avail] <- coef(fit_c)[avail]
    vc_out <- matrix(NA_real_, length(regs), length(regs),
                     dimnames = list(regs, regs))
    vc_out[avail, avail] <- vc
    list(coef = cf_out, vcov = vc_out, iter = k)
  }, error = function(e) { cat("    GNLS error:", e$message, "\n"); NULL })
}
gnls <- run_gnls(fit_ppml, panel, fml_3d, REGRESSORS)
if (!is.null(gnls)) {
  cf <- gnls$coef; se <- sqrt(diag(gnls$vcov)); pv <- 2 * pnorm(-abs(cf / se))
  results[["GNLS"]] <- tibble(estimator = "GNLS", variable = names(cf),
                              coef = round(cf, 6), se = round(se, 6),
                              pval = round(pv, 6), stars = add_stars(pv))
}

# ---- [8] OLS-3D: log-linearised comparator --------------------------------
cat("  [8/9] OLS-3D...\n")
fit_ols <- tryCatch(
  feols(fml_ols, data = panel, vcov = ~pair_id),
  error = function(e) { cat("    OLS error:", e$message, "\n"); NULL }
)
results[["OLS_3D"]] <- get_results(fit_ols, "OLS-3D")

# ---- [9] 2W-PPML: two-way comparator --------------------------------------
cat("  [9/9] 2W-PPML...\n")
fit_2w <- tryCatch(
  fepois(fml_2w, data = panel, vcov = ~pair_id, glm.iter = 200, nthreads = 1),
  error = function(e) { cat("    2W error:", e$message, "\n"); NULL }
)
results[["PPML_2W"]] <- get_results(fit_2w, "2W-PPML")

# ---- Assemble the wide Table 3 (full-panel block) -------------------------
results_all <- bind_rows(results)

t3 <- results_all %>%
  mutate(cell = mapply(fmt_cell, coef, se, stars)) %>%
  select(estimator, variable, cell) %>%
  pivot_wider(names_from = estimator, values_from = cell) %>%
  mutate(Variable = dplyr::case_when(
    variable == "rta"         ~ "RTA",
    variable == "ln_dist"     ~ "ln(distance)",
    variable == "contig"      ~ "Contiguity",
    variable == "comlang_off" ~ "Common language",
    variable == "comcol"      ~ "Common coloniser",
    TRUE ~ variable)) %>%
  select(Variable, any_of(EST_ORDER))

t3 <- bind_rows(
  t3,
  tibble(Variable = "Observations",     SPJ = format(nrow(panel), big.mark = ",")),
  tibble(Variable = "Zero trade share",
         SPJ = sprintf("%.1f%%", 100 * mean(panel$trade_kusd == 0))),
  tibble(Variable = "Variance-to-mean ratio (non-zero)",
         SPJ = format(round(var(nz) / mean(nz)), big.mark = ","))
)

write.csv(results_all, "output/table3_full_panel_long.csv", row.names = FALSE)
write.csv(t3,          "output/table3_full_panel_wide.csv", row.names = FALSE)

cat("\n  RTA coefficient across estimators:\n")
print(results_all %>% filter(variable == "rta") %>%
        select(estimator, coef, se, stars), n = 20)

cat("\n03_main_estimation.R complete: wrote Table 3 full-panel block.\n\n")
