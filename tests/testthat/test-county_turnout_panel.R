# Structural / integrity checks for the bundled county_turnout_panel dataset.
# These assert the invariants that data-raw/30_validate.R checks at build time,
# so a corrupted or accidentally-regenerated .rda is caught by R CMD check.

data("county_turnout_panel", package = "pintervals")
p <- county_turnout_panel

test_that("schema matches the documented column set and order", {
  expect_s3_class(p, "tbl_df")
  expect_identical(
    names(p),
    c("state", "county", "fips", "election_year",
      "turnout", "total_votes", "cvap", "votes_dem", "votes_rep", "votes_other",
      "dem_share", "rep_share", "other_share",
      "turnout_county_mean", "turnout_electorate",
      "dem_share_county_mean", "dem_share_electorate",
      "rep_share_county_mean", "rep_share_electorate",
      "other_share_county_mean", "other_share_electorate",
      "total_population", "nonwhite_pct", "foreignborn_pct", "female_pct",
      "age29andunder_pct", "age65andolder_pct", "median_hh_inc",
      "clf_unemploy_pct", "lesscollege_pct", "lesshs_pct", "rural_pct",
      "ruralurban_cc", "division", "region", "geo_group",
      "longitude", "latitude")
  )
  expect_type(p$fips, "character")
  expect_true(all(grepl("^[0-9]{5}$", p$fips)))
  expect_true(is.factor(p$division) && is.factor(p$region) && is.factor(p$geo_group))
})

test_that("the panel key is unique (no duplicate county-year rows)", {
  expect_equal(anyDuplicated(p[c("fips", "election_year")]), 0L)
})

test_that("exactly the five presidential elections are present", {
  expect_setequal(unique(p$election_year), c(2008L, 2012L, 2016L, 2020L, 2024L))
})

test_that("county counts per year match the investigated values", {
  tab <- as.integer(table(p$election_year))
  # native geography per wave is 3114 / 3114 / 3113 / 3113 / 3114 (2016 loses
  # Bedford city VA -> Bedford County; 2024 gains 1 from CT 8 counties -> 9
  # planning regions); then 6 / 3 / 4 / 6 / 6 turnout > 1 micro-county-years are
  # dropped. See the codebook / R/data.R.
  expect_identical(tab, c(3108L, 3111L, 3109L, 3107L, 3108L))
})

test_that("turnout: values > 1 are dropped, NA is retained", {
  # micro-county county-years whose ballots exceed the noisy 5-year ACS CVAP
  # estimate are removed in data-raw/20_assemble_panel.R.
  expect_true(all(p$turnout[!is.na(p$turnout)] <= 1))
  expect_true(all(p$turnout[!is.na(p$turnout)] > 0.05))
  expect_true(anyNA(p$turnout))
})

test_that("turnout NA is limited to the documented cases", {
  na_rows <- p[is.na(p$turnout), ]
  expect_lte(nrow(na_rows), 14L)
  expect_true(all(na_rows$county == "Kalawao" | na_rows$state == "Connecticut"))
})

test_that("turnout is reproducible from total_votes and cvap", {
  expect_type(p$total_votes, "integer")
  expect_type(p$cvap, "integer")
  ok <- !is.na(p$turnout)
  expect_equal(round(p$total_votes[ok] / p$cvap[ok], 6), p$turnout[ok])
  expect_identical(which(is.na(p$total_votes)), which(is.na(p$turnout)))
  expect_false(anyNA(p$cvap))
  expect_true(all(p$cvap > 0))
})

test_that("votes by party reconcile to total_votes and shares sum to ~1", {
  for (cc in c("votes_dem", "votes_rep", "votes_other")) {
    expect_type(p[[cc]], "integer")
    expect_true(all(p[[cc]][!is.na(p[[cc]])] >= 0), info = cc)
    expect_identical(which(is.na(p[[cc]])), which(is.na(p$turnout)), info = cc)
  }
  ok <- !is.na(p$turnout)
  resid <- (p$votes_dem + p$votes_rep + p$votes_other - p$total_votes)[ok]
  # exact for 2008-2020; <= 13 tiny 2024 counties differ, none by much
  expect_lt(max(abs(resid)), 3000)
  expect_lte(sum(abs(resid) > 50), 15L)
  for (cc in c("dem_share", "rep_share", "other_share")) {
    v <- p[[cc]][!is.na(p[[cc]])]
    expect_true(all(v >= 0 & v <= 1), info = cc)
  }
  ssum <- (p$dem_share + p$rep_share + p$other_share)[ok]
  expect_true(all(abs(ssum - 1) < 0.02))
})

test_that("per-election aggregates are constant within a wave and plausible", {
  agg_cols <- c("turnout_county_mean", "turnout_electorate",
                "dem_share_county_mean", "dem_share_electorate",
                "rep_share_county_mean", "rep_share_electorate",
                "other_share_county_mean", "other_share_electorate")
  for (cc in agg_cols) {
    expect_false(anyNA(p[[cc]]), info = cc)
    per_wave <- tapply(p[[cc]], p$election_year, function(x) length(unique(x)))
    expect_true(all(per_wave == 1), info = cc)
  }
  agg <- unique(p[c("election_year", agg_cols)])
  expect_true(all(agg$turnout_electorate > 0.5 & agg$turnout_electorate < 0.75))
  expect_true(all(abs(agg$dem_share_electorate + agg$rep_share_electorate +
                      agg$other_share_electorate - 1) < 0.005))
})

test_that("percentage covariates lie in [0, 100]", {
  for (cc in c("nonwhite_pct", "foreignborn_pct", "female_pct", "age29andunder_pct",
               "age65andolder_pct", "clf_unemploy_pct", "lesscollege_pct",
               "lesshs_pct", "rural_pct")) {
    v <- p[[cc]][!is.na(p[[cc]])]
    expect_true(all(v >= 0 & v <= 100), info = cc)
  }
})

test_that("identifiers, codes and coordinates are in range", {
  expect_true(all(p$total_population > 0))
  expect_true(all(p$median_hh_inc[!is.na(p$median_hh_inc)] > 0))
  expect_true(all(p$ruralurban_cc[!is.na(p$ruralurban_cc)] %in% 1:9))
  expect_false(anyNA(p$ruralurban_cc))
  expect_true(all(p$longitude > -180 & p$longitude < -60))
  expect_true(all(p$latitude > 15 & p$latitude < 72))
  expect_false(anyNA(p$division) || anyNA(p$region) || anyNA(p$geo_group))
})

test_that("covariate missingness stays within documented bounds", {
  expect_lte(sum(is.na(p$rural_pct)), 9L)      # CT 2024 planning regions
  expect_lte(sum(is.na(p$median_hh_inc)), 5L)  # Census-suppressed micro-counties
})

test_that("Alaska and the island territories are excluded", {
  expect_false(any(substr(p$fips, 1, 2) %in% c("02", "60", "66", "69", "72", "78")))
})
