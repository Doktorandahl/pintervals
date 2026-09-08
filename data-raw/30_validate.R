# ---------------------------------------------------------------------------
# 30_validate.R -- Post-assembly validation report for county_turnout_panel
# ---------------------------------------------------------------------------
# Prints a human-readable report and stops on hard failures. County-years with
# turnout > 1 are dropped in 20_assemble_panel.R; this script asserts none
# remain. NA outcomes are retained and reported, never dropped. The same
# invariants are asserted against the shipped .rda in
# tests/testthat/test-county_turnout_panel.R.
# ---------------------------------------------------------------------------

if (!exists("county_turnout_panel")) source(file.path("data-raw", "20_assemble_panel.R"))

p    <- county_turnout_panel
ok   <- function(...) cat("  [ok]   ", sprintf(...), "\n")
flag <- function(...) cat("  [flag] ", sprintf(...), "\n")
hard <- function(cond, msg) if (!isTRUE(cond)) stop("VALIDATION FAILURE: ", msg, call. = FALSE)

cat("\n=== county_turnout_panel validation ===\n\n")

## 1. shape / schema -------------------------------------------------------------
expected_cols <- c("state", "county", "fips", "election_year", "turnout",
  "total_population", "nonwhite_pct", "foreignborn_pct", "female_pct",
  "age29andunder_pct", "age65andolder_pct", "median_hh_inc", "clf_unemploy_pct",
  "lesscollege_pct", "lesshs_pct", "rural_pct", "ruralurban_cc",
  "division", "region", "geo_group", "longitude", "latitude")
hard(identical(names(p), expected_cols), "column set / order does not match the spec")
ok("schema: %d columns in expected order", ncol(p))
hard(is.character(p$fips) && all(grepl("^[0-9]{5}$", p$fips)), "fips must be 5-char strings")
ok("fips: all 5-character strings")

## 2. no duplicate county-year rows -------------------------------------------
dups <- p[duplicated(p[c("fips", "election_year")]) |
          duplicated(p[c("fips", "election_year")], fromLast = TRUE), ]
hard(nrow(dups) == 0, sprintf("%d duplicate fips-year rows", nrow(dups)))
ok("uniqueness: 0 duplicate (fips, election_year) rows")

## 3. county counts per year -- investigate, don't assert constant -----------
cnt <- as.data.frame(table(p$election_year))
names(cnt) <- c("election_year", "n_counties")
cat("\n  county rows per election year:\n")
for (i in seq_len(nrow(cnt))) cat(sprintf("    %s : %d\n", cnt[i, 1], cnt[i, 2]))
delta <- diff(cnt$n_counties)
cat("\n  wave-to-wave change: ", paste(sprintf("%+d", delta), collapse = "  "), "\n")
# Drivers of the differences (all documented in the codebook):
#  * native geography per wave: 3114 / 3114 / 3113 / 3113 / 3114 county-equivalents
#    (2016 loses Bedford (ind.) city VA -> Bedford County; 2024 gains 1 from the
#     Connecticut 8 legacy counties -> 9 planning regions change)
#  * then the turnout > 1 micro-county-years dropped in 20_assemble_panel.R:
#    6 / 3 / 4 / 6 / 6, leaving 3108 / 3111 / 3109 / 3107 / 3108
expected_counts <- c(`2008` = 3108, `2012` = 3111, `2016` = 3109, `2020` = 3107, `2024` = 3108)
hard(identical(as.integer(cnt$n_counties), unname(as.integer(expected_counts[as.character(cnt$election_year)]))),
     "per-year county counts differ from the investigated/expected values")
ok("county counts: match the investigated expectations (see codebook)")

## 4. turnout in [0, 1] -- > 1 dropped upstream, NA retained ------------------
n_na  <- sum(is.na(p$turnout))
oob   <- p[!is.na(p$turnout) & (p$turnout < 0 | p$turnout > 1), ]
cat("\n  turnout: ", n_na, " NA | ", nrow(oob), " outside [0, 1]\n")
if (n_na) {
  na_tbl <- as.data.frame(p[is.na(p$turnout), c("election_year", "state", "county")])
  flag("turnout NA for %d rows (Kalawao Co. HI, all waves; Connecticut planning regions, 2024):", n_na)
  print(na_tbl, row.names = FALSE)
}
# turnout > 1 is dropped in 20_assemble_panel.R -- assert none slipped through
hard(nrow(oob) == 0, sprintf("%d rows with turnout outside [0, 1] survived the drop", nrow(oob)))
hard(all(p$turnout[!is.na(p$turnout)] > 0.05), "turnout has implausibly low values")
hard(n_na <= 14, "more NA turnout than the documented 14")
ok("turnout: <= 1 for every non-NA row; NA only for the documented cases")

## 5. covariate plausibility ------------------------------------------------
pct_cols <- c("nonwhite_pct", "foreignborn_pct", "female_pct", "age29andunder_pct",
              "age65andolder_pct", "clf_unemploy_pct", "lesscollege_pct",
              "lesshs_pct", "rural_pct")
for (cc in pct_cols) {
  v <- p[[cc]]; v <- v[!is.na(v)]
  hard(all(v >= 0 & v <= 100), sprintf("%s outside [0, 100]", cc))
}
ok("percentage covariates: all within [0, 100]")
hard(all(p$total_population > 0), "non-positive total_population")
hard(all(p$median_hh_inc[!is.na(p$median_hh_inc)] > 0), "non-positive median_hh_inc")
hard(sum(is.na(p$median_hh_inc)) <= 5,
     "more than 5 NA median_hh_inc (expected: Census-suppressed in a couple of micro-counties)")
hard(all(p$ruralurban_cc[!is.na(p$ruralurban_cc)] %in% 1:9), "ruralurban_cc not in 1:9")
hard(all(p$longitude > -180 & p$longitude < -60), "longitude out of CONUS+HI range")
hard(all(p$latitude > 15 & p$latitude < 72), "latitude out of range")
ok("population, income, RUCC codes, coordinates: all in range")

## 6. bounded missingness in covariates ------------------------------------
covar_na <- vapply(p[setdiff(expected_cols, c("turnout"))], function(x) sum(is.na(x)), integer(1))
covar_na <- covar_na[covar_na > 0]
if (length(covar_na)) {
  cat("\n  covariate NA counts:\n")
  for (nm in names(covar_na)) flag("%-16s %d NA", nm, covar_na[[nm]])
}
hard(sum(is.na(p$rural_pct)) <= 9, "unexpected rural_pct NA (> the 9 CT-2024 rows)")
hard(sum(is.na(p$ruralurban_cc)) == 0, "unexpected ruralurban_cc NA")
ok("covariate missingness: within documented bounds")

## 7. geo classification completeness -------------------------------------
hard(!anyNA(p$division) && !anyNA(p$region) && !anyNA(p$geo_group),
     "missing division / region / geo_group")
ok("division / region / geo_group: complete for every row")

## 8. object footprint ----------------------------------------------------
sz <- as.numeric(object.size(p))
cat(sprintf("\n  in-memory size: %.2f MB\n", sz / 1024^2))
tf <- tempfile(fileext = ".rda")
save(p, file = tf, compress = "xz")
ser_kb <- file.size(tf) / 1024
unlink(tf)
cat(sprintf("  serialized (xz): %.0f KB\n", ser_kb))
hard(is.finite(ser_kb) && ser_kb < 1.5 * 1024, "serialized panel larger than 1.5 MB")
ok("footprint: comfortably within CRAN data-size expectations")

cat("\n=== validation complete: no hard failures ===\n")
invisible(TRUE)
