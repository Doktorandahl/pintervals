# ---------------------------------------------------------------------------
# 20_assemble_panel.R -- Join the tidy sources into county_turnout_panel
# ---------------------------------------------------------------------------
# Depends on objects from 10_transform.R (+ 07_fetch_kc_split.R). Produces:
#   county_turnout_panel   the bundled tibble
#   panel_build_notes      a list of everything reconciled / flagged
#
# Geography policy: the ACS 5-year county list for each wave is the spine
# ("native geography per wave"). Alaska and the island territories are
# excluded (Alaska has no county/borough-level presidential returns; this
# also matches the original `county_turnout`). Where another source uses a
# different FIPS vintage than the spine, it is remapped; residual join
# mismatches are left as NA and recorded in `panel_build_notes`, never silently
# dropped. County-years whose computed turnout exceeds 1 (recorded ballots above
# the noisy 5-year ACS CVAP estimate, a handful of micro-counties) ARE dropped
# from the shipped panel and recorded in `panel_build_notes$turnout_over_1_dropped`.
# ---------------------------------------------------------------------------

if (!exists("acs_src"))  source(file.path("data-raw", "10_transform.R"))
if (!exists("kc_split")) source(file.path("data-raw", "07_fetch_kc_split.R"))
library(stringr)

notes <- list()

## ==========================================================================
## 1. Spine: ACS county list per wave, 50 states + DC
## ==========================================================================
territory_st <- c("60", "66", "69", "72", "78")            # AS GU MP PR VI
ak_st        <- "02"

spine <- acs_src |>
  mutate(st = substr(fips, 1, 2)) |>
  filter(!st %in% c(territory_st, ak_st)) |>
  transmute(
    fips, election_year,
    county = str_trim(str_remove(county_api_name, ",[^,]+$")),
    state  = str_trim(str_extract(county_api_name, "[^,]+$"))
  ) |>
  mutate(county = str_trim(str_remove(county, "\\s+(County|Parish)$")))

notes$spine_rows_per_year <- spine |> count(election_year, name = "n_counties")
notes$excluded_geographies <- paste(
  "Alaska (state 02): no county/borough-level presidential returns in MEDSL;",
  "American Samoa (60), Guam (66), N. Mariana Is. (69), Puerto Rico (72),",
  "U.S. Virgin Islands (78): not states, no presidential electors / not in MEDSL."
)

## ==========================================================================
## 2. Turnout numerator: MEDSL votes + Kansas City redistribution
## ==========================================================================
kc_fips <- c("29037", "29047", "29095", "29165")           # Cass, Clay, Jackson, Platte

votes <- turnout_src |>
  transmute(fips = fips_raw, election_year, total_votes)

# Kansas City rows: FIPS 2938000 (<=2020) or 3600000 / mis-coded 36000 (2024).
kc_rows <- turnout_src |>
  filter(medsl_county == "KANSAS CITY") |>
  transmute(election_year, kc_votes = total_votes)
notes$kansas_city <- turnout_src |>
  filter(medsl_county == "KANSAS CITY") |>
  transmute(election_year, fips_in_medsl = fips_raw, kc_votes = total_votes)

kc_add <- purrr::map_dfr(seq_len(nrow(vintages)), function(i) {
  ey <- vintages$election_year[i]
  w  <- kc_split[[as.character(vintages$dec_vintage[i])]] |>
        filter(fips %in% kc_fips)
  kv <- kc_rows$kc_votes[kc_rows$election_year == ey]
  if (!length(kv)) return(tibble())
  tibble(fips = w$fips, election_year = ey, kc_extra_votes = w$kc_weight * kv)
})

votes <- votes |>
  filter(!fips %in% c("2938000", "3600000", "36000")) |>   # drop the KC unit itself
  left_join(kc_add, by = c("fips", "election_year")) |>
  mutate(total_votes = total_votes + dplyr::coalesce(kc_extra_votes, 0)) |>
  select(fips, election_year, total_votes)

# FIPS remap of the numerator onto the spine vintage
#   Shannon -> Oglala Lakota (46113 -> 46102) for waves whose spine uses 46102
#   Bedford (ind.) city -> Bedford County (51515 -> 51019) for the 2016 wave
votes <- votes |>
  mutate(fips = dplyr::case_when(
    fips == "46113" & election_year %in% c(2016, 2020, 2024) ~ "46102",
    fips == "51515" & election_year == 2016                  ~ "51019",
    TRUE ~ fips
  )) |>
  group_by(fips, election_year) |>
  summarise(total_votes = sum(total_votes), .groups = "drop")

## ==========================================================================
## 3. Denominator: CVAP (already on the spine's ACS vintage)
## ==========================================================================
den <- cvap_src |> select(fips, election_year, cvap)

## ==========================================================================
## 4. Rural / urban continuum code (USDA ERS)
## ==========================================================================
# rename Shannon; drop the standalone Bedford-city RUCC row (2016 wave keeps
# the Bedford *County* code). Connecticut 2020: the 2016-2020 ACS spine uses
# the 8 legacy counties but RUCC 2023 is tabulated on the 9 planning regions,
# so CT-2020 takes its code from RUCC 2013 (legacy counties) instead.
rucc_2013_ct <- read_rucc(2013L) |>
  filter(substr(fips, 1, 2) == "09") |>
  mutate(election_year = 2020L)

rucc <- rucc_src |>
  mutate(fips = dplyr::if_else(fips == "46113" & election_year %in% c(2016, 2020, 2024),
                               "46102", fips)) |>
  filter(!(fips == "51515" & election_year == 2016)) |>   # 2016 spine has only 51019
  filter(!(substr(fips, 1, 2) == "09" & election_year == 2020)) |>   # replace CT 2020
  bind_rows(rucc_2013_ct) |>
  distinct(fips, election_year, .keep_all = TRUE) |>
  select(fips, election_year, ruralurban_cc)

notes$rucc_vintage_by_wave <- vintages |>
  transmute(election_year, rucc_vintage,
            note = dplyr::case_when(
              election_year %in% c(2012, 2016) ~ "RUCC 2013 shared by 2012 & 2016 waves (no ERS revision between)",
              election_year %in% c(2020, 2024) ~ "RUCC 2023 shared by 2020 & 2024 waves (no ERS revision between)",
              TRUE ~ "RUCC 2003"))
notes$rucc_connecticut_2020 <- "CT 2020 wave uses RUCC 2013 (legacy 8 counties); RUCC 2023 switched to the 9 planning regions and does not join the 2016-2020 ACS geography."

## ==========================================================================
## 5. rural_pct from the nearest decennial census (P2 urban/rural)
## ==========================================================================
# Count-level FIPS reconciliation, then recompute the percentage:
#   Shannon -> Oglala Lakota rename; Bedford city + county combined (2016 wave).
rural <- rural_src |>
  mutate(fips = dplyr::case_when(
    fips == "46113" & election_year %in% c(2016, 2020, 2024) ~ "46102",
    fips == "51515" & election_year == 2016                  ~ "51019",
    TRUE ~ fips
  )) |>
  group_by(fips, election_year) |>
  summarise(pop_total_dec = sum(pop_total_dec),
            pop_rural_dec = sum(pop_rural_dec), .groups = "drop") |>
  mutate(rural_pct = 100 * pop_rural_dec / pop_total_dec) |>
  select(fips, election_year, rural_pct)

notes$rural_pct_source <- paste(
  "rural_pct is held fixed from the nearest decennial census (2010 Census P2 for the",
  "2008/2012/2016 waves; 2020 Census P2 for the 2020/2024 waves) because the ACS",
  "carries no urban/rural classification. Connecticut 2024: the 2020 Census P2 is on",
  "the 8 legacy counties while the 2020-2024 ACS spine is on the 9 planning regions,",
  "so CT-2024 rural_pct is NA.")

## ==========================================================================
## 6. County centroid longitude / latitude (Census Gazetteer)
## ==========================================================================
geo <- geo_src |> select(fips, election_year, longitude, latitude)

## ==========================================================================
## 7. division / region / geo_group  (inherited state-level classification)
## ==========================================================================
data_env <- new.env(); load("data/county_turnout.rda", envir = data_env)
geo_class <- data_env$county_turnout |>
  dplyr::distinct(state, division, region, geo_group) |>
  as_tibble()
stopifnot(setequal(setdiff(unique(spine$state), "District of Columbia"),
                   setdiff(geo_class$state, "District of Columbia")))
notes$geo_class_source <- paste(
  "division and region are the U.S. Census Bureau's regions and divisions;",
  "geo_group is a pintervals-specific coarse grouping. Both are a fixed",
  "state-level classification carried over verbatim from the original",
  "`county_turnout` dataset.")

## ==========================================================================
## 8. Assemble
## ==========================================================================
county_turnout_panel <- spine |>
  left_join(votes, by = c("fips", "election_year")) |>
  left_join(den,   by = c("fips", "election_year")) |>
  left_join(select(acs_src, -county_api_name), by = c("fips", "election_year")) |>
  left_join(rural, by = c("fips", "election_year")) |>
  left_join(rucc,  by = c("fips", "election_year")) |>
  left_join(geo,   by = c("fips", "election_year")) |>
  left_join(geo_class, by = "state") |>
  mutate(turnout = total_votes / cvap) |>
  transmute(
    state, county, fips,
    election_year = as.integer(election_year),
    turnout          = round(turnout, 6),
    total_population = as.integer(round(total_population)),
    nonwhite_pct      = round(nonwhite_pct, 4),
    foreignborn_pct   = round(foreignborn_pct, 4),
    female_pct        = round(female_pct, 4),
    age29andunder_pct = round(age29andunder_pct, 4),
    age65andolder_pct = round(age65andolder_pct, 4),
    median_hh_inc     = as.integer(round(median_hh_inc)),
    clf_unemploy_pct  = round(clf_unemploy_pct, 4),
    lesscollege_pct   = round(lesscollege_pct, 4),
    lesshs_pct        = round(lesshs_pct, 4),
    rural_pct         = round(rural_pct, 4),
    ruralurban_cc     = as.integer(ruralurban_cc),
    division  = factor(division),
    region    = factor(region),
    geo_group = factor(geo_group),
    longitude = round(longitude, 5),
    latitude  = round(latitude, 5)
  ) |>
  arrange(election_year, fips)

## ==========================================================================
## 9. Drop implausible turnout (> 1); flag rows the joins could not complete
## ==========================================================================
# turnout > 1 shows up only in micro-counties where recorded presidential
# ballots exceed the noisy 5-year ACS CVAP estimate. These county-years are
# removed from the shipped panel (rows with NA turnout are retained). The
# dropped rows are recorded here, never discarded silently.
notes$turnout_over_1_dropped <- county_turnout_panel |>
  filter(!is.na(turnout), turnout > 1) |>
  select(election_year, fips, state, county, turnout) |>
  arrange(dplyr::desc(turnout))

county_turnout_panel <- county_turnout_panel |>
  filter(is.na(turnout) | turnout <= 1)

notes$missing_turnout <- county_turnout_panel |>
  filter(is.na(turnout)) |>
  select(election_year, fips, state, county)
notes$turnout_out_of_unit_interval <- county_turnout_panel |>
  filter(!is.na(turnout), turnout < 0 | turnout > 1) |>
  select(election_year, fips, state, county, turnout)
notes$missing_any_covariate <- county_turnout_panel |>
  filter(dplyr::if_any(c(rural_pct, ruralurban_cc, longitude, median_hh_inc), is.na)) |>
  count(election_year, name = "n_rows_with_a_missing_covariate")

panel_build_notes <- notes

message(sprintf("assembled county_turnout_panel: %d rows, %d cols, years %s",
                nrow(county_turnout_panel), ncol(county_turnout_panel),
                paste(sort(unique(county_turnout_panel$election_year)), collapse = ", ")))
message(sprintf("  NA turnout: %d (retained) | turnout > 1: %d (dropped)",
                nrow(notes$missing_turnout), nrow(notes$turnout_over_1_dropped)))
