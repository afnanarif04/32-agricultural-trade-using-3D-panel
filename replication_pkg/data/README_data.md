# Data documentation

**File:** `BiasCorr3DPPML_Agri_MergedPanel.csv`
**Rows:** 40,560 (directed exporter–importer–year cells)
**Panel:** 40 exporters × 40 importers × 26 years (1995–2020), self-pairs excluded.

The panel merges the CEPII structural gravity backbone, BACI HS92 bilateral
agricultural trade flows, and World Bank macro indicators into one
analysis-ready file. All variables below are already processed; no further
cleaning is required to reproduce the paper.

---

## Variables

### Identifiers and dimensions

| Column | Description | Unit / type |
|--------|-------------|-------------|
| `iso3_o` | Exporter (origin) ISO3 code | string |
| `iso3_d` | Importer (destination) ISO3 code | string |
| `year` | Year | 1995–2020 |
| `pair` | Exporter–importer pair label | string |
| `pair_id` | Pair fixed-effect index (mu_ij) | integer |
| `exp_time_id` | Exporter-time fixed-effect index (lambda_it) | integer |
| `imp_time_id` | Importer-time fixed-effect index (xi_jt) | integer |
| `trend` | Linear time trend | integer |

### Outcome

| Column | Description | Unit |
|--------|-------------|------|
| `y_ijt` | Bilateral agricultural trade value | USD thousands |
| `ln_y_ijt` | log(1 + `y_ijt`) | — |

### Gravity regressors and dyadic controls

| Column | Description | Unit / type |
|--------|-------------|-------------|
| `fta_wto` | Regional trade agreement / WTO indicator (the "RTA" variable) | 0/1 |
| `dist`, `distcap` | Bilateral distance (population-weighted; capital-to-capital) | km |
| `ln_dist` | log distance | — |
| `contig` | Shared land border | 0/1 |
| `comlang_off` | Common official language | 0/1 |
| `comcol` | Common coloniser post-1945 | 0/1 |
| `col45` | Colonial relationship post-1945 | 0/1 |

### Macro covariates (exporter `_o`, importer `_d`)

| Column | Description | Unit |
|--------|-------------|------|
| `gdp_fin_o`, `gdp_fin_d` | GDP (gravity dataset) | current USD |
| `gdpcap_fin_o`, `gdpcap_fin_d` | GDP per capita (gravity dataset) | current USD |
| `pop_fin_o`, `pop_fin_d` | Population (gravity dataset) | persons |
| `ln_gdp_o`, `ln_gdp_d`, `ln_gdpcap_o`, `ln_gdpcap_d`, `ln_pop_o`, `ln_pop_d` | Logs of the above | — |
| `GDP_current_USD_o/d`, `GDPpc_current_USD_o/d`, `Population_o/d`, `Trade_openness_pct_GDP_o/d` | World Bank WDI macro indicators | USD / persons / % of GDP |

### Membership and event flags

| Column | Description | Unit |
|--------|-------------|------|
| `wto_o`, `wto_d` | WTO membership (exporter / importer) | 0/1 |
| `eu_o`, `eu_d` | EU membership (exporter / importer) | 0/1 |
| `crisis_2008` | 2008 global financial crisis flag | 0/1 |
| `covid_2020` | 2020 COVID flag | 0/1 |
| `post_paris` | Post-2015 Paris Agreement flag | 0/1 |

---

## Data sources

| Source | Variables | URL | Accessed |
|--------|-----------|-----|----------|
| CEPII Gravity Dataset | distance, contiguity, language, colonial ties, gravity GDP/pop | http://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=8 | 2026-06-29 |
| CEPII BACI (HS92, V202601) | bilateral agricultural trade flows `y_ijt` | http://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37 | 2026-06-29 |
| World Bank World Development Indicators | GDP, GDP per capita, population, trade openness | https://databank.worldbank.org/source/world-development-indicators | 2026-06-29 |

Agricultural trade is defined as BACI HS92 chapters 01–24 aggregated to the
directed dyad-year. Country coverage is the 40 economies listed in the paper's
appendix; self-pairs (i = j) are excluded.

---

## Sample

- **N (exporters) = N (importers) = 40**, **T = 26** (1995–2020).
- **Observations:** 40,560 directed dyad-year cells; 39,570 enter estimation
  after `fepois` removes singleton fixed-effect groups.
- **Frequency:** annual.
- **Zero-trade share:** 2.8% of directed dyad-year cells.
- **Variance-to-mean ratio (non-zero trade):** on the order of several million,
  indicating severe overdispersion (motivating the negative-binomial DGP in the
  Monte Carlo).

---

## Notes

- **Missing values.** The macro covariates contain occasional gaps; they are
  not used as regressors in the reported specifications (which rely on the
  exporter-time and importer-time fixed effects to absorb all monadic variation)
  and therefore do not affect replication of the tables.
- **Panel balance.** The dyad-year grid is balanced by construction; zeros are
  genuine zero-trade observations, not missing data.
- **Cross-section-invariant variables.** `dist`, `ln_dist`, `contig`,
  `comlang_off`, `comcol`, and `col45` are time-invariant within a pair and are
  therefore absorbed by the pair fixed effect in the three-way specifications;
  they are identified only in the two-way comparator.
