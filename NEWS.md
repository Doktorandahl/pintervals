
# pintervals 1.3.1

* `county_turnout_panel` now excludes the 25 micro-county-years whose recorded
  presidential ballots exceeded the noisy 5-year ACS CVAP estimate (turnout > 1),
  leaving 15,543 rows. Rows with an unknown (`NA`) turnout are still retained.
  The filtering is performed by the build pipeline in `data-raw/`.

# pintervals 1.3.0

* New dataset `county_turnout_panel`: a tidy five-wave county-level panel of
  U.S. presidential-election turnout and ACS demographic context for the 2008,
  2012, 2016, 2020 and 2024 elections (15,568 rows). Turnout is votes cast
  divided by Citizen Voting-Age Population (CVAP). All demographic covariates
  are re-derived directly from the U.S. Census Bureau API (not IPUMS NHGIS, to
  comply with NHGIS's non-redistribution terms) and match the vintage of each
  election. Ships with a full codebook at
  `system.file("extdata", "county_turnout_panel_codebook.md", package = "pintervals")`
  and a reproducible build pipeline in `data-raw/`.
* `county_turnout` (the original single-election 2016 cross-section) is
  **superseded** by `county_turnout_panel` but **retained** for backward
  compatibility; it is still used by the package vignette and examples. Its data
  are unchanged.
* Documentation fix: `?county_turnout` previously stated that `turnout` was
  "total votes cast divided by total population". It is in fact votes divided by
  the Citizen Voting-Age Population (CVAP); the help page has been corrected.
  `county_turnout_panel$turnout` uses this same definition and reproduces the
  2016 `county_turnout` values for 3,103 of the 3,107 shared counties. The Details
  section of `?county_turnout` (the `predicted_turnout` reproduction code) also
  now renders — it was previously dropped because of an unescaped `%` in the
  Rd source.

# pintervals 1.2.0

* Updated Readme
* Updated Vignette

# pintervals 1.1.4

* Bugfixes and minor improvements.
