
## Update

This is an update of pintervals from version 1.1.1 (current CRAN release) to
1.3.3. Main changes since 1.1.1:

* New dataset `county_turnout_panel`: a five-wave county-level panel of U.S.
  presidential-election turnout, vote shares and ACS demographics (2008–2024).
  All data are public domain or CC0 (U.S. Census Bureau API, MEDSL election
  returns), with a codebook in `inst/extdata` and the build pipeline in
  `data-raw/`. The existing `county_turnout` dataset is retained unchanged
  for backward compatibility.
* Vignette rewritten to match the manuscript under review at the Journal of
  Statistical Software; `ranger` added to Suggests.
* Documentation fixes, `importFrom("stats", "lm.fit")` added to NAMESPACE,
  and bug fixes.

## Maintainer email change

The maintainer email address has changed from david.randahl@pcr.uu.se to david.randahl@fhs.se


## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New maintainer: the maintainer email address has changed (see above).

## Reverse dependencies

There are currently no reverse dependencies.
