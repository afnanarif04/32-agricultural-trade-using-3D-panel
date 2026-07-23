# Replication Package

Replication code and data for a study of finite-sample bias correction in
three-way fixed-effects Poisson gravity estimation, applied to bilateral
agricultural commodity trade. The package reproduces every table in the paper:
a Monte Carlo evaluation of the proposed correction and an empirical
application to a 40-country agricultural trade panel.

*(Author, title, and affiliation are omitted for anonymous review and will be
added on acceptance.)*

---

## Software requirements

- **R** version 4.1.0 or later (tested on 4.3.x).
- The following packages (installed automatically by `01_packages_data.R` if
  missing):

```r
install.packages(c("data.table", "dplyr", "tidyr", "fixest"))
```

`fixest` (>= 0.11) provides the `fepois`/`feols` estimators used throughout.
`RcppParallel`, if present, is used only to pin execution to a single thread
for reproducibility.

---

## Repository structure

```
.
├── README.md
├── LICENSE
├── data/
│   ├── README_data.md
│   └── BiasCorr3DPPML_Agri_MergedPanel.csv   ← processed analysis panel
├── code/
│   ├── 00_run_all.R          master runner
│   ├── 01_packages_data.R    packages, data load, shared objects/helpers
│   ├── 02_descriptives.R     summary statistics, fixed-effect counts
│   ├── 03_main_estimation.R  nine estimators, full panel  (Table 3, full)
│   ├── 04_robustness_spj.R   SPJ across four subsamples    (Table 3, robust.)
│   └── 05_monte_carlo.R      81-cell simulation            (Table 1)
└── output/
    └── .gitkeep              (result CSVs are written here)
```

---

## Data

- **Panel:** 40 exporters × 40 importers × 26 years (1995–2020).
- **Observations:** 40,560 directed dyad-year cells (39,570 after singleton
  removal inside `fepois`).
- **Unit of observation:** directed exporter–importer–year.
- Full variable and source documentation is in `data/README_data.md`.

---

## How to replicate

**Option 1 — run everything.** Set the working directory to the repository
root (the folder containing `code/`, `data/`, `output/`) and source the master
runner:

```r
setwd("/path/to/repository-root")
source("code/00_run_all.R")
```

**Option 2 — run one table at a time.** Each script sources
`01_packages_data.R` itself, so any can be run alone:

```r
setwd("/path/to/repository-root")
source("code/03_main_estimation.R")   # Table 3, full-panel block
source("code/04_robustness_spj.R")    # Table 3, robustness block
source("code/05_monte_carlo.R")       # Table 1
```

All result files are written to `output/`.

---

## Output → paper table map

| Script | Output file(s) | Paper table |
|--------|----------------|-------------|
| `02_descriptives.R` | `table_descriptives_summary.csv`, `table_descriptives_fe_counts.csv` | Section 4.1 in-text figures |
| `03_main_estimation.R` | `table3_full_panel_wide.csv`, `table3_full_panel_long.csv` | Table 3 (full-panel block) |
| `04_robustness_spj.R` | `table3_robustness_wide.csv`, `table3_robustness_long.csv` | Table 3 (robustness block) |
| `05_monte_carlo.R` | `summary_bias_rmse.csv`, `summary_size.csv`, `summary_power_spj.csv` | Table 1 |

---

## Estimators

| Code label | Estimator | Reference |
|------------|-----------|-----------|
| `SPJ` | Split-panel-jackknife bias correction (proposed) | Dhaene & Jochmans (2015) framework |
| `3D-PPML` | Uncorrected three-way Poisson PMLE | Weidner & Zylkin (2021) |
| `WZ-Analytical` | Leave-one-pair-out approximation to the WZ correction | Weidner & Zylkin (2021) |
| `Polyad` | Polyad estimator (placeholder; public code unavailable) | Resende et al. (2025) |
| `Bootstrap` | Pair-bootstrap bias correction | Higgins & Jochmans (2024) |
| `Heckman-2S` | Selection two-step | Mnasri & Nechi (2021) |
| `GNLS` | Iterated generalised nonlinear least squares | Mnasri & Nechi (2021) |
| `OLS-3D` | Log-linearised OLS with three-way fixed effects | — |
| `2W-PPML` | Two-way Poisson comparator (pair + year) | — |

Five estimators (`SPJ`, `3D-PPML`, `WZ-Analytical`, `GNLS`, `OLS-3D`) are
reported as columns in Table 3; the remainder are discussed in Section 4.5.

---

## Notes

- **Seed.** A global seed of 42 is set in `00_run_all.R` and each script; the
  Monte Carlo sets a distinct seed per grid cell (`1000 + cell index`).
- **Runtime.** The empirical scripts (`02`–`04`) run in a few minutes on a
  standard laptop. `03_main_estimation.R` is dominated by the `B = 200`
  bootstrap; reduce `B_BOOT` at the top of that script for a faster check.
  `05_monte_carlo.R` (81 cells × 2,000 replications) is the slowest component
  and can take several hours single-threaded; set `QUICK <- TRUE` at the top of
  that script for a small-grid structural check.
- **Threads.** Execution is pinned to a single thread for reproducibility. To
  speed up on a multi-core machine, raise `nthreads` in the `fepois`/`feols`
  calls, at the cost of exact numerical reproducibility.
- **Tested on** R 4.3.x, Linux and Windows.
