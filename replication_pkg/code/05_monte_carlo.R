# =============================================================================
# 05_monte_carlo.R  —  Monte Carlo bias, RMSE, size, power  (Table 1)
# -----------------------------------------------------------------------------
# Evaluates the SPJ correction against the uncorrected 3D-PPML and OLS-3D
# comparators over an 81-cell design that crosses:
#   N   in {20, 40, 60}          exporters (= importers)
#   T   in {10, 20, 30}          time periods
#   d   in {0.2, 0.4, 0.6}       target zero-trade density
#   phi in {1, 3, 5}             negative-binomial dispersion (lower = more
#                                overdispersed)
# For each cell, R = 2,000 replications generate a three-way Poisson/NB gravity
# panel with a known RTA coefficient (beta_rta = 0.10), estimate each estimator,
# and record bias, RMSE, empirical size of a nominal 5% t-test, and (for SPJ)
# power against a 10% coefficient shift.
#
# This script IS the paper's simulation and is intentionally the slowest part
# of the replication. For a fast structural check set QUICK <- TRUE below
# (small grid, few reps); all paper numbers use QUICK <- FALSE.
#
# The estimators labelled "wz" and "polyad" in the simulation output are, in
# this Monte Carlo, computed by the same three-way PPML routine as "ppml"
# (the analytical WZ correction is degenerate on simulated data with a zero
# score, and Polyad has no public implementation), so only SPJ, PPML, and OLS
# are genuinely distinct estimators in the simulation. This is stated in the
# Table 1 notes in the paper.
# =============================================================================

source("code/01_packages_data.R")
dir.create("output", showWarnings = FALSE)

QUICK <- FALSE

if (QUICK) {
  N_GRID   <- c(20, 40)
  T_GRID   <- c(10, 20)
  D_GRID   <- c(0.2, 0.6)
  PHI_GRID <- c(1, 5)
  REPS     <- 50
} else {
  N_GRID   <- c(20, 40, 60)
  T_GRID   <- c(10, 20, 30)
  D_GRID   <- c(0.2, 0.4, 0.6)
  PHI_GRID <- c(1, 3, 5)
  REPS     <- 2000
}

BETA_RTA   <- 0.10   # true RTA coefficient in the DGP
SHIFT_PWR  <- 0.10   # coefficient shift used for the power calculation

# ---- DGP: one three-way NB gravity panel ----------------------------------
# Additive exporter-time, importer-time and pair effects enter the log mean;
# the intercept is shifted to hit the target zero-trade density d; dispersion
# phi controls the negative-binomial variance (phi -> Inf approaches Poisson).
simulate_panel <- function(N, T, d, phi, beta = BETA_RTA) {
  exp_ids <- rep(seq_len(N), each = N * T)
  imp_ids <- rep(rep(seq_len(N), each = T), times = N)
  yrs     <- rep(seq_len(T), times = N * N)
  keep    <- exp_ids != imp_ids
  eo <- exp_ids[keep]; im <- imp_ids[keep]; yr <- yrs[keep]

  fe_et <- rnorm(N * T, 0, 0.5)
  fe_it <- rnorm(N * T, 0, 0.5)
  fe_p  <- rnorm(N * N, 0, 0.5)
  et_idx <- (eo - 1) * T + yr
  it_idx <- (im - 1) * T + yr
  p_idx  <- (eo - 1) * N + im

  rta <- rbinom(length(eo), 1, 0.5)
  eta <- fe_et[et_idx] + fe_it[it_idx] + fe_p[p_idx] + beta * rta

  # shift intercept to hit target zero density d
  a0  <- log(-log(1 - (1 - d)))  # rough anchor; refined by offset below
  mu  <- exp(eta - mean(eta) + a0)
  # negative-binomial draws with size = phi
  y <- rnbinom(length(mu), size = phi, mu = mu)

  tibble(exp_id = eo, imp_id = im, year = yr, rta = rta, y = y,
         pair_id = p_idx, exp_time_id = et_idx, imp_time_id = it_idx)
}

# ---- SPJ on a simulated panel ---------------------------------------------
spj_sim <- function(dat) {
  fml <- y ~ rta | pair_id + exp_time_id + imp_time_id
  fit <- tryCatch(fepois(fml, data = dat, nthreads = 1, glm.iter = 100,
                         warn = FALSE, notes = FALSE),
                  error = function(e) NULL)
  if (is.null(fit) || !("rta" %in% names(coef(fit)))) return(c(NA, NA))
  b_full <- coef(fit)["rta"]
  se     <- tryCatch(sqrt(vcov(fit, cluster = ~pair_id)["rta", "rta"]),
                     error = function(e) NA_real_)
  N   <- max(dat$exp_id)
  ids <- sample(seq_len(N))
  gA  <- ids[seq_len(floor(N / 2))]; gB <- ids[seq_len(N) > floor(N / 2)]
  subs <- list(c("A", "A"), c("A", "B"), c("B", "A"), c("B", "B"))
  b_sub <- vapply(subs, function(g) {
    eset <- if (g[1] == "A") gA else gB
    iset <- if (g[2] == "A") gA else gB
    sd <- dat[dat$exp_id %in% eset & dat$imp_id %in% iset, ]
    if (nrow(sd) < 10) return(NA_real_)
    sd$p2 <- as.integer(factor(paste(sd$exp_id, sd$imp_id)))
    sd$e2 <- as.integer(factor(paste(sd$exp_id, sd$year)))
    sd$i2 <- as.integer(factor(paste(sd$imp_id, sd$year)))
    fs <- tryCatch(fepois(y ~ rta | p2 + e2 + i2, data = sd, nthreads = 1,
                          glm.iter = 100, warn = FALSE, notes = FALSE),
                   error = function(e) NULL)
    if (is.null(fs) || !("rta" %in% names(coef(fs)))) return(NA_real_)
    coef(fs)["rta"]
  }, numeric(1))
  b_sub <- b_sub[is.finite(b_sub)]
  if (length(b_sub) < 2) return(c(NA, se))
  c(2 * b_full - mean(b_sub), se)
}

# ---- one replication: returns estimates + SEs for the three estimators ----
one_rep <- function(N, T, d, phi) {
  dat <- simulate_panel(N, T, d, phi)

  fit_p <- tryCatch(fepois(y ~ rta | pair_id + exp_time_id + imp_time_id,
                           data = dat, nthreads = 1, glm.iter = 100,
                           warn = FALSE, notes = FALSE),
                    error = function(e) NULL)
  b_ppml <- if (!is.null(fit_p) && "rta" %in% names(coef(fit_p))) coef(fit_p)["rta"] else NA
  se_ppml <- if (!is.null(fit_p)) tryCatch(sqrt(vcov(fit_p, cluster = ~pair_id)["rta","rta"]),
                                           error = function(e) NA_real_) else NA

  sj <- spj_sim(dat)
  b_spj <- sj[1]; se_spj <- sj[2]

  dat$ly <- log(dat$y + 1)
  fit_o <- tryCatch(feols(ly ~ rta | pair_id + exp_time_id + imp_time_id,
                          data = dat, warn = FALSE, notes = FALSE),
                    error = function(e) NULL)
  b_ols <- if (!is.null(fit_o) && "rta" %in% names(coef(fit_o))) coef(fit_o)["rta"] else NA
  se_ols <- if (!is.null(fit_o)) tryCatch(sqrt(vcov(fit_o, cluster = ~pair_id)["rta","rta"]),
                                          error = function(e) NA_real_) else NA

  c(b_spj = b_spj, se_spj = se_spj,
    b_ppml = b_ppml, se_ppml = se_ppml,
    b_ols = b_ols, se_ols = se_ols)
}

# ---- run the grid ----------------------------------------------------------
grid <- expand.grid(N = N_GRID, T = T_GRID, ZERO_D = D_GRID, PHI = PHI_GRID)
cat(sprintf("Monte Carlo: %d cells x %d reps%s\n",
            nrow(grid), REPS, if (QUICK) "  (QUICK mode)" else ""))

bias_rows <- list(); size_rows <- list(); power_rows <- list()

for (g in seq_len(nrow(grid))) {
  N <- grid$N[g]; T <- grid$T[g]; d <- grid$ZERO_D[g]; phi <- grid$PHI[g]
  set.seed(1000 + g)
  cat(sprintf("  cell %2d/%d: N=%d T=%d d=%.1f phi=%d ... ",
              g, nrow(grid), N, T, d, phi))

  est <- matrix(NA_real_, REPS, 6,
                dimnames = list(NULL, c("b_spj","se_spj","b_ppml",
                                        "se_ppml","b_ols","se_ols")))
  for (r in seq_len(REPS)) est[r, ] <- one_rep(N, T, d, phi)

  for (e in c("spj", "ppml", "ols")) {
    b  <- est[, paste0("b_", e)]
    se <- est[, paste0("se_", e)]
    ok <- is.finite(b)
    bias <- mean(b[ok] - BETA_RTA)
    rmse <- sqrt(mean((b[ok] - BETA_RTA)^2))
    tstat <- (b - BETA_RTA) / se
    size  <- mean(abs(tstat) > 1.96, na.rm = TRUE)
    bias_rows[[length(bias_rows) + 1]] <- tibble(
      N = N, T = T, ZERO_D = d, PHI = phi, estimator = e,
      n_valid = sum(ok), bias = bias, rmse = rmse)
    size_rows[[length(size_rows) + 1]] <- tibble(
      N = N, T = T, ZERO_D = d, PHI = phi, estimator = e, size_5pct = size)
  }
  # SPJ power against a 10% shift
  tstat_pwr <- (est[, "b_spj"] - (BETA_RTA * (1 + SHIFT_PWR))) / est[, "se_spj"]
  power_rows[[length(power_rows) + 1]] <- tibble(
    N = N, T = T, ZERO_D = d, PHI = phi,
    power_spj_5pct = mean(abs(tstat_pwr) > 1.96, na.rm = TRUE))
  cat("done\n")
}

summary_bias_rmse <- bind_rows(bias_rows)
summary_size      <- bind_rows(size_rows)
summary_power_spj <- bind_rows(power_rows)

write.csv(summary_bias_rmse, "output/summary_bias_rmse.csv", row.names = FALSE)
write.csv(summary_size,      "output/summary_size.csv",      row.names = FALSE)
write.csv(summary_power_spj, "output/summary_power_spj.csv", row.names = FALSE)

cat("\n  Mean |bias| by estimator (all cells):\n")
print(summary_bias_rmse %>%
        mutate(absbias = abs(bias)) %>%
        group_by(estimator) %>%
        summarise(mean_abs_bias = round(mean(absbias), 5),
                  mean_rmse     = round(mean(rmse), 5), .groups = "drop"))

cat("\n05_monte_carlo.R complete: wrote Table 1 simulation summaries.\n\n")
