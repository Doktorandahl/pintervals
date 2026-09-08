# data-raw/ — building the bundled datasets

## `county_turnout_panel`

A 5-wave (2008, 2012, 2016, 2020, 2024) county-level panel of U.S. presidential
turnout and ACS demographic context. Rebuilds `data/county_turnout_panel.rda`.

### How to regenerate from scratch

1. Get a free Census Bureau API key: <https://api.census.gov/data/key_signup.html>
2. Add it to `~/.Renviron` (then restart R):

   ```
   CENSUS_API_KEY=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   ```

3. From the package root:

   ```sh
   Rscript data-raw/county_turnout_panel.R
   ```

The run is **idempotent**: every raw download is cached under
`data-raw/cache/` (git-ignored, ~260 MB) and reused on later runs. Delete a
cache sub-directory to force that source to re-download. Fetching is fully
separated from transformation:

| script | stage | needs key | network |
|:-------|:------|:---------:|:--------|
| `00_config.R` | election-year → vintage map; cache + API helpers | – | – |
| `01_fetch_medsl.R` | MEDSL county presidential returns (Dataverse API) | no | yes |
| `02_fetch_cvap.R` | CVAP special tabulations (census.gov/rdo) | no | yes |
| `03_fetch_acs.R` | ACS 5-year detailed tables (Census API) | **yes** | yes |
| `04_fetch_rucc.R` | USDA ERS Rural-Urban Continuum Codes | no | yes |
| `05_fetch_geo.R` | Census Gazetteer county centroids | no | yes |
| `06_fetch_decennial.R` | decennial census urban/rural split (Census API) | **yes** | yes |
| `07_fetch_kc_split.R` | Kansas City place-by-county population weights (Census API) | **yes** | yes |
| `10_transform.R` | each raw source → one tidy county-year frame | – | – |
| `20_assemble_panel.R` | joins + FIPS/geography harmonisation → the panel | – | – |
| `30_validate.R` | schema / duplicates / turnout range / county counts | – | – |
| `county_turnout_panel.R` | orchestrator: runs all of the above, writes the `.rda` | **yes** | yes |

Supporting files: `fips_crosswalk.csv` (documented geography changes with
sources), `acs_variables_resolved.csv` (written by `03`; every ACS table pulled,
per vintage, for traceability).

Runtime is dominated by downloads (~260 MB) and ~70 Census API calls; a cold run
takes a few minutes, a warm run (cache populated) about 15 seconds.

### Exact upstream citations

* **MEDSL vote returns** —
  MIT Election Data and Science Lab, 2024,
  *County Presidential Election Returns 2000–2024*,
  <https://doi.org/10.7910/DVN/VOQCHQ>, Harvard Dataverse, V20. **License: CC0 1.0.**
  Downloaded via the Dataverse data-access API. The dataset content is CC0, but
  Harvard Dataverse attaches a *required* download guestbook (name / email /
  institution / position); `01_fetch_medsl.R` satisfies it with a guestbook
  response POST (no account or token needed — edit `gb_response` in that script
  to record your own details).

* **CVAP (turnout denominator)** —
  U.S. Census Bureau, *Citizen Voting-Age Population (CVAP) Special Tabulation
  From the American Community Survey*, 2006–2010, 2008–2012, 2012–2016,
  2016–2020 and 2020–2024 5-year files,
  <https://www.census.gov/programs-surveys/decennial-census/about/voting-rights/cvap.html>.
  U.S. Government work — public domain.

* **ACS demographic covariates** —
  U.S. Census Bureau, *American Community Survey 5-Year Estimates*, Detailed
  Tables (B01001, B01003, B03002, B05002, B15002, B15003, B19013, B23001,
  B23025), obtained directly from the Census Bureau API,
  <https://www.census.gov/data/developers/data-sets/acs-5year.html>.
  U.S. Government work — public domain. **Not** sourced via IPUMS NHGIS (see
  below).

* **rural_pct** —
  U.S. Census Bureau, 2010 Census Summary File 1 table P2 and 2020 Census
  Demographic and Housing Characteristics File table P2 ("Urban and Rural"), via
  the Census Bureau API. U.S. Government work — public domain.

* **Rural-Urban Continuum Codes** —
  U.S. Department of Agriculture, Economic Research Service, *Rural-Urban
  Continuum Codes*, 2003, 2013 and 2023 releases,
  <https://www.ers.usda.gov/data-products/rural-urban-continuum-codes/>.
  U.S. Government work — public domain.

* **County centroids** —
  U.S. Census Bureau, *Gazetteer Files* (county records, 2012 / 2013 / 2016 /
  2020 / 2024),
  <https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.html>.
  U.S. Government work — public domain.

* **County boundary-change documentation** (for `fips_crosswalk.csv`) —
  U.S. Census Bureau, *Substantial Changes to Counties and County Equivalent
  Entities: 1970–Present*,
  <https://www.census.gov/programs-surveys/geography/technical-documentation/county-changes.html>;
  and *Change to County Equivalents in the State of Connecticut*, 87 FR 34235
  (2022-06-06).

### Licensing / data-consistency notes discovered during assembly

* **NHGIS avoidance is load-bearing.** IPUMS NHGIS prohibits redistribution of
  extracted data. All ACS/decennial covariates therefore come straight from the
  Census Bureau API. Validated: the reconstruction reproduces the NHGIS-sourced
  `county_turnout` to floating-point tolerance.
* **MEDSL Harvard Dataverse guestbook** — not a licence term, but the build must
  POST a guestbook response to download. Handled in `01_fetch_medsl.R`.
* **`county_turnout` help page was wrong about its own `turnout`.** It claims
  votes / total population; the values are actually votes / CVAP (verified: the
  panel's independently-computed votes/CVAP matches `county_turnout$turnout` to
  < 0.001 for 3,103 / 3,107 counties). The help page is corrected in this
  release; the `county_turnout` data are untouched.
* **MEDSL revised some 2016 county returns** (Gila/Graham/Greenlee AZ, Keya Paha
  NE) between the version used to build `county_turnout` and v20 — the panel
  uses the corrected returns, so its 2016 turnout differs from `county_turnout`
  for those four counties.
* **MEDSL 2024 Kansas City FIPS is mis-coded** as `36000` (name still
  "KANSAS CITY") where 2000–2020 use `2938000`. `20_assemble_panel.R` detects it
  by name and redistributes it correctly; flagged here as an upstream error.
* **Alaska** has no county/borough-level presidential returns in MEDSL (reported
  by state-house district). Alaska is excluded from the panel.
* **Connecticut** switched from 8 counties to 9 planning regions across the
  window, and every source adopted the change on a different schedule. Kept
  native per wave; consequences documented in the codebook and the roxygen.
* **Micro-county turnout > 1** — 25 county-years where actual ballots exceed the
  5-year ACS CVAP point estimate. Retained and flagged, per the build request.

## `county_turnout` and `elections`

The original `county_turnout` (single-election 2016 cross-section) and
`elections` datasets predate this directory and have no build script here.
`county_turnout` is superseded by `county_turnout_panel` but retained for
backward compatibility (it is still used by the package vignette).
