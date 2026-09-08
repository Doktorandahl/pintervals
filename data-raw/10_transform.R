# ---------------------------------------------------------------------------
# 10_transform.R -- Turn each cached raw source into a tidy county-year frame
# ---------------------------------------------------------------------------
# Consumes only files under data-raw/cache/ (populated by 01-06_fetch_*.R).
# Produces, in the calling environment:
#   turnout_src   fips_raw, election_year, total_votes,
#                 votes_dem, votes_rep, votes_other                (MEDSL)
#   cvap_src      fips, election_year, cvap                       (CVAP special tab)
#   acs_src       fips, election_year, <12 covariates>            (ACS 5-year)
#   rucc_src      fips, election_year, ruralurban_cc              (USDA ERS)
#   rural_src     fips, election_year, rural_pct                  (decennial P2)
#   geo_src       fips, election_year, longitude, latitude        (Census Gazetteer)
#   acs_spine     fips, election_year, county_api_name            (ACS county list)
#
# All geographic reconciliation (FIPS crosswalking, Kansas City, Connecticut
# 2024, Alaska) happens later, in 20_assemble_panel.R.
# ---------------------------------------------------------------------------

if (!exists("vintages")) source(file.path("data-raw", "00_config.R"))
library(stringr)

target_years <- vintages$election_year

## ==========================================================================
## 1. MEDSL -> turnout numerator + votes by party (Dem / Rep / other)
## ==========================================================================
medsl_raw <- readr::read_tsv(cache_path("medsl", "countypres_2000-2024.tab"),
                             show_col_types = FALSE, guess_max = 1e6)

medsl_pres <- medsl_raw |> filter(year %in% target_years, office == "US PRESIDENT")

# total_votes: the county-wide presidential ballot count, repeated on every
# candidate/mode row; it is unique within county-year (verified in build).
tot_src <- medsl_pres |>
  group_by(election_year = as.integer(year), state, county_name, county_fips) |>
  summarise(total_votes = dplyr::first(totalvotes), .groups = "drop")

# votes by party. Complications handled here:
#   * MEDSL only breaks 2008-2016 into DEMOCRAT / REPUBLICAN / OTHER; 2020 adds
#     GREEN + LIBERTARIAN, 2024 adds LIBERTARIAN. Everything that is not
#     DEMOCRAT/REPUBLICAN is folded into a single "other" group so the panel is
#     consistent across all five waves.
#   * some 2020/2024 counties report by voting mode with no "TOTAL" row -> sum
#     the mode rows; where a "TOTAL" row exists it already is the sum -> use it.
#   * 2024 carries junk `party == NA` rows ("TOTAL VOTES CAST", "UNDERVOTES",
#     "OVERVOTES", "SPOILED") that would double-count -> dropped.
party_by_cand <- medsl_pres |>
  filter(!is.na(party), !is.na(county_fips)) |>
  mutate(mode = dplyr::coalesce(mode, "TOTAL"),
         candidatevotes = dplyr::coalesce(candidatevotes, 0)) |>
  group_by(election_year = as.integer(year), county_fips, party, candidate) |>
  summarise(votes = if (any(mode == "TOTAL")) sum(candidatevotes[mode == "TOTAL"])
                    else sum(candidatevotes),
            .groups = "drop")

party_wide <- party_by_cand |>
  mutate(pg = dplyr::case_when(party == "DEMOCRAT"   ~ "dem",
                               party == "REPUBLICAN" ~ "rep",
                               TRUE                  ~ "other")) |>
  group_by(election_year, county_fips) |>
  summarise(votes_dem   = sum(votes[pg == "dem"]),
            votes_rep   = sum(votes[pg == "rep"]),
            votes_other = sum(votes[pg == "other"]),
            .groups = "drop")

turnout_src <- tot_src |>
  left_join(party_wide, by = c("election_year", "county_fips")) |>
  mutate(votes_dem   = dplyr::coalesce(votes_dem, 0),
         votes_rep   = dplyr::coalesce(votes_rep, 0),
         votes_other = dplyr::coalesce(votes_other, 0),
         fips_raw = dplyr::case_when(
           is.na(county_fips)     ~ NA_character_,
           county_fips == 2938000 ~ "2938000",               # Kansas City, MO (own unit)
           TRUE                   ~ str_pad(county_fips, 5, "left", "0")
         ))

# rows deliberately not carried into the county panel (documented in codebook):
#   * county_fips == NA : state-level write-in / UOCAVA / "Federal Precinct" buckets
#   * state == ALASKA   : MEDSL reports AK by state-house district, not borough
turnout_drop_note <- turnout_src |>
  filter(is.na(fips_raw) | state == "ALASKA") |>
  count(state, reason = dplyr::if_else(state == "ALASKA",
                                       "no borough-level returns (reported by house district)",
                                       "state-level ballot bucket, county_fips missing"),
        name = "rows")

turnout_src <- turnout_src |>
  filter(!is.na(fips_raw), state != "ALASKA") |>
  select(fips_raw, election_year, medsl_state = state, medsl_county = county_name,
         total_votes, votes_dem, votes_rep, votes_other)

## ==========================================================================
## 2. CVAP -> turnout denominator
## ==========================================================================
read_cvap_county <- function(label) {
  f <- list.files(cache_path("cvap", label), pattern = "(?i)county\\.csv$",
                  full.names = TRUE, recursive = TRUE)[1]
  d <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
  names(d) <- tolower(names(d))
  d |>
    filter(.data$lntitle == "Total") |>                 # all race/ethnicity groups
    transmute(fips = str_extract(.data$geoid, "\\d{5}$"),
              cvap = as.numeric(.data$cvap_est))
}

cvap_src <- purrr::map_dfr(seq_len(nrow(vintages)), function(i) {
  read_cvap_county(vintages$cvap_label[i]) |>
    mutate(election_year = vintages$election_year[i], .after = fips)
})

## ==========================================================================
## 3. ACS -> demographic covariates (+ county spine)
## ==========================================================================
# All ACS variables used here are non-negative counts or a median dollar
# amount. The API encodes "not computable / suppressed" as large negative jam
# values (-666666666, -999999999, ...) and occasionally -1; map every negative
# to NA so it propagates rather than corrupting a ratio.
acs_num <- function(df, v) {
  x <- suppressWarnings(as.numeric(df[[v]]))
  x[x < 0] <- NA_real_
  x
}

# Sum a set of estimate columns row-wise, tolerating a single-column selection.
row_sum_vars <- function(df, vars) {
  if (length(vars) == 0) stop("row_sum_vars: no variables selected", call. = FALSE)
  m <- vapply(vars, function(v) acs_num(df, v), numeric(nrow(df)))
  if (is.null(dim(m))) return(m)          # length-1 selection -> vector
  rowSums(m)
}

# Normalise an ACS variable label to a stable path form.
# 2009-2018 vintages: "Estimate!!Total!!Male!!Under 5 years"
# 2019+   vintages:   "Estimate!!Total:!!Male:!!Under 5 years"
# -> both become      "Total!!Male!!Under 5 years"
norm_lbl <- function(x) {
  x |>
    stringr::str_remove("^\\s*Estimate!!") |>
    stringr::str_replace_all(":!!", "!!") |>
    stringr::str_remove(":$") |>
    stringr::str_trim()
}

lbl_ages <- function(norm_label) {
  if (grepl("and over", norm_label)) return(c(min = suppressWarnings(min(as.numeric(stringr::str_extract_all(norm_label, "\\d+")[[1]]))), max = Inf))
  n <- as.numeric(stringr::str_extract_all(norm_label, "\\d+")[[1]])
  if (!length(n)) return(c(min = NA_real_, max = NA_real_))
  c(min = min(n), max = max(n))
}

transform_acs_wave <- function(end_year) {
  ey  <- as.character(end_year)
  rd  <- function(tb) readRDS(cache_path("acs", ey, paste0(tb, ".rds")))
  lb  <- function(tb) {
    d <- readRDS(cache_path("acs", ey, paste0(tb, "_labels.rds")))
    d$norm <- norm_lbl(d$label)
    d
  }

  b01003 <- rd("B01003"); b03002 <- rd("B03002"); b05002 <- rd("B05002")
  b19013 <- rd("B19013"); b01001 <- rd("B01001")
  edu_tb <- if (end_year <= 2011L) "B15002" else "B15003"
  emp_tb <- if (end_year <= 2011L) "B23001" else "B23025"
  edu <- rd(edu_tb); emp <- rd(emp_tb)
  edu_l <- lb(edu_tb); emp_l <- lb(emp_tb); age_l <- lb("B01001")

  # every table for a vintage is pulled with the same geo predicate, so row
  # order is identical; assert it before using positional alignment.
  key <- function(d) fips5(d$state, d$county)
  stopifnot(identical(key(b01003), key(b01001)), identical(key(b03002), key(b01001)),
            identical(key(b05002), key(b01001)), identical(key(b19013), key(b01001)),
            identical(key(edu),    key(b01001)), identical(key(emp),    key(b01001)))

  # --- age shares from B01001 (label-driven; structure identical all vintages)
  age_detail <- age_l |>
    filter(stringr::str_count(norm, "!!") == 2) |>        # "Total!!<sex>!!<age>" leaves
    mutate(amin = vapply(norm, function(z) lbl_ages(z)[["min"]], numeric(1)),
           amax = vapply(norm, function(z) lbl_ages(z)[["max"]], numeric(1)))
  v_u29 <- age_detail$name[age_detail$amax <= 29]
  v_65p <- age_detail$name[age_detail$amin >= 65]
  stopifnot(length(v_u29) == 18, length(v_65p) == 12)     # 9 male + 9 female / 6 + 6
  pop_total <- acs_num(b01001, "B01001_001E")
  age29andunder_pct <- 100 * row_sum_vars(b01001, v_u29) / pop_total
  age65andolder_pct <- 100 * row_sum_vars(b01001, v_65p) / pop_total
  female_pct        <- 100 * acs_num(b01001, "B01001_026E") / pop_total

  # --- education from B15003 (>=2012) or B15002 (2010) -----------------------
  below_hs <- c("No schooling completed", "Nursery school", "Kindergarten",
                "1st grade", "2nd grade", "3rd grade", "4th grade", "5th grade",
                "6th grade", "7th grade", "8th grade", "9th grade", "10th grade",
                "11th grade", "12th grade, no diploma",
                # B15002 (2010) collapses some of the above:
                "Nursery to 4th grade", "5th and 6th grade", "7th and 8th grade")
  college  <- c("Bachelor's degree", "Master's degree",
                "Professional school degree", "Doctorate degree")
  if (edu_tb == "B15003") {
    edu_tot  <- acs_num(edu, "B15003_001E")
    lesshs_v <- edu_l$name[edu_l$norm %in% paste0("Total!!", below_hs)]
    coll_v   <- edu_l$name[edu_l$norm %in% paste0("Total!!", college)]
    stopifnot(length(lesshs_v) == 15, length(coll_v) == 4)
  } else {                                                # B15002: sex-split
    edu_tot  <- acs_num(edu, "B15002_001E")
    lesshs_v <- edu_l$name[edu_l$norm %in% as.vector(outer(c("Total!!Male!!", "Total!!Female!!"), below_hs, paste0))]
    coll_v   <- edu_l$name[edu_l$norm %in% as.vector(outer(c("Total!!Male!!", "Total!!Female!!"), college, paste0))]
    stopifnot(length(lesshs_v) == 16, length(coll_v) == 8)
  }
  lesshs_pct      <- 100 * row_sum_vars(edu, lesshs_v) / edu_tot
  lesscollege_pct <- 100 * (1 - row_sum_vars(edu, coll_v) / edu_tot)

  # --- unemployment from B23025 (>=2012) or B23001 (2010) -------------------
  if (emp_tb == "B23025") {
    clf_unemploy_pct <- 100 * acs_num(emp, "B23025_005E") / acs_num(emp, "B23025_003E")
  } else {
    un_v <- emp_l$name[stringr::str_ends(emp_l$norm, "!!Unemployed")]
    em_v <- emp_l$name[stringr::str_ends(emp_l$norm, "!!Employed")]
    un   <- row_sum_vars(emp, un_v)
    em   <- row_sum_vars(emp, em_v)
    clf_unemploy_pct <- 100 * un / (un + em)              # civilian labour force = emp + unemp
  }

  tibble(
    fips              = fips5(b01001$state, b01001$county),
    election_year     = vintages$election_year[match(end_year, vintages$acs_end)],
    county_api_name   = b01003$NAME,
    total_population  = acs_num(b01003, "B01003_001E"),
    nonwhite_pct      = 100 * (1 - acs_num(b03002, "B03002_003E") / acs_num(b03002, "B03002_001E")),
    foreignborn_pct   = 100 * acs_num(b05002, "B05002_013E") / acs_num(b05002, "B05002_001E"),
    female_pct        = female_pct,
    age29andunder_pct = age29andunder_pct,
    age65andolder_pct = age65andolder_pct,
    median_hh_inc     = acs_num(b19013, "B19013_001E"),
    clf_unemploy_pct  = clf_unemploy_pct,
    lesscollege_pct   = lesscollege_pct,
    lesshs_pct        = lesshs_pct
  )
}

acs_src <- purrr::map_dfr(vintages$acs_end, transform_acs_wave)

acs_spine <- acs_src |> select(fips, election_year, county_api_name)

## ==========================================================================
## 4. USDA ERS Rural-Urban Continuum Codes
## ==========================================================================
read_rucc <- function(vintage) {
  if (vintage == 2023L) {
    d <- readr::read_csv(cache_path("rucc", "rucc_2023.csv"), show_col_types = FALSE)
    d |> filter(.data$Attribute == "RUCC_2023") |>
      transmute(fips = str_pad(.data$FIPS, 5, "left", "0"),
                ruralurban_cc = as.numeric(.data$Value))
  } else if (vintage == 2013L) {
    d <- readxl::read_excel(cache_path("rucc", "rucc_2013.xls"))
    transmute(d, fips = str_pad(.data$FIPS, 5, "left", "0"),
              ruralurban_cc = as.numeric(.data$RUCC_2013))
  } else {
    d <- readxl::read_excel(cache_path("rucc", "rucc_2003.xls"))
    transmute(d, fips = str_pad(.data$`FIPS Code`, 5, "left", "0"),
              ruralurban_cc = as.numeric(.data$`2003 Rural-urban Continuum Code`))
  }
}

rucc_src <- purrr::map_dfr(seq_len(nrow(vintages)), function(i) {
  read_rucc(vintages$rucc_vintage[i]) |>
    mutate(election_year = vintages$election_year[i], .after = fips)
})

## ==========================================================================
## 5. Decennial rural / urban split -> rural_pct
## ==========================================================================
rural_by_year <- list(
  `2010` = readRDS(cache_path("decennial", "rural_2010.rds")),
  `2020` = readRDS(cache_path("decennial", "rural_2020.rds"))
)
rural_src <- purrr::map_dfr(seq_len(nrow(vintages)), function(i) {
  cy <- as.character(vintages$dec_vintage[i])
  rural_by_year[[cy]] |>
    transmute(fips, election_year = vintages$election_year[i],
              pop_total_dec, pop_rural_dec, rural_pct)   # counts kept for FIPS merges
})

## ==========================================================================
## 6. Census Gazetteer -> county centroid longitude / latitude
## ==========================================================================
gaz_year <- c(`2008` = 2012L, `2012` = 2013L, `2016` = 2016L, `2020` = 2020L, `2024` = 2024L)
read_gaz <- function(gy) {
  d <- readr::read_tsv(cache_path("geo", sprintf("%d_Gaz_counties_national.txt", gy)),
                       show_col_types = FALSE, name_repair = "minimal")
  names(d) <- trimws(names(d))
  transmute(d,
            fips = str_pad(.data$GEOID, 5, "left", "0"),
            longitude = as.numeric(.data$INTPTLONG),
            latitude  = as.numeric(.data$INTPTLAT))
}
geo_src <- purrr::map_dfr(names(gaz_year), function(ey) {
  read_gaz(gaz_year[[ey]]) |> mutate(election_year = as.integer(ey), .after = fips)
})

message("10_transform: ",
        nrow(turnout_src), " turnout | ", nrow(cvap_src), " cvap | ",
        nrow(acs_src), " acs | ", nrow(rucc_src), " rucc | ",
        nrow(rural_src), " rural | ", nrow(geo_src), " geo rows")
