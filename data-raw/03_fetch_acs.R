# ---------------------------------------------------------------------------
# 03_fetch_acs.R -- Demographic covariates from ACS 5-year detailed tables
# ---------------------------------------------------------------------------
# Source : U.S. Census Bureau, American Community Survey 5-year estimates,
#          DETAILED TABLES, pulled DIRECTLY from the Census Bureau API.
#          U.S. Government work, public domain.
#
#          Deliberately NOT sourced via IPUMS NHGIS: NHGIS's terms of use
#          prohibit redistributing extracted data outside narrow journal
#          exceptions, so NHGIS-derived data cannot be bundled in a CRAN
#          package. The identical ACS tables from the Census Bureau's own API
#          carry no such restriction.
#
# Access : Census Bureau API (needs CENSUS_API_KEY in ~/.Renviron).
# Output : data-raw/cache/acs/<endYear>/<TABLE>.rds        (raw wide counts)
#          data-raw/cache/acs/<endYear>/<TABLE>_labels.rds (variable metadata)
#          data-raw/acs_variables_resolved.csv             (traceability index)
#
# ACS window per election wave (window whose END year is nearest the election):
#   2008 -> 2006-2010 (vintage 2010) | 2012 -> 2008-2012 | 2016 -> 2012-2016
#   2020 -> 2016-2020                | 2024 -> 2020-2024
#
# Table availability note: B15003 (education) and B23025 (employment status)
# do not exist in the 2010 vintage; that wave uses B15002 and B23001 instead.
# Every other table is consistent across all five vintages.
# ---------------------------------------------------------------------------

if (!exists("vintages")) source(file.path("data-raw", "00_config.R"))

# tables needed in EVERY vintage
acs_tables_common <- c(
  total_population = "B01003",   # total population
  race             = "B03002",   # Hispanic or Latino origin by race
  nativity         = "B05002",   # place of birth by citizenship status
  income           = "B19013",   # median household income (12-mo, infl-adj)
  agesex           = "B01001"    # sex by age
)
# education + employment: table id depends on vintage
acs_table_education  <- function(end) if (end <= 2011L) "B15002" else "B15003"
acs_table_employment <- function(end) if (end <= 2011L) "B23001" else "B23025"

acs_raw     <- list()   # acs_raw[[endYear]][[table]]  -> wide tibble
acs_labels  <- list()   # acs_labels[[endYear]][[table]] -> name/label/concept
resolved    <- list()   # rows for acs_variables_resolved.csv

for (i in seq_len(nrow(vintages))) {
  end   <- vintages$acs_end[i]
  ey    <- as.character(end)
  path  <- sprintf("%d/acs/acs5", end)
  outd  <- cache_path("acs", ey)
  dir.create(outd, showWarnings = FALSE, recursive = TRUE)

  tables <- c(acs_tables_common,
              education  = acs_table_education(end),
              employment = acs_table_employment(end))

  message(sprintf("\n  ACS %d (%d-%d): %s",
                  end, vintages$acs_start[i], end, paste(tables, collapse = ", ")))

  acs_raw[[ey]]    <- list()
  acs_labels[[ey]] <- list()

  for (tb in tables) {
    lab_cache <- file.path(outd, sprintf("%s_labels.rds", tb))
    dat_cache <- file.path(outd, sprintf("%s.rds", tb))

    labs <- census_group_vars(path, tb, cache = lab_cache)
    # estimate columns only (names end in E, exclude margins of error "M")
    est_vars <- labs$name[grepl("E$", labs$name)]

    # NAME gives the geography's official label for this exact ACS vintage
    # (e.g. "Autauga County, Alabama" or, for CT 2024, a planning-region name).
    dat <- census_get(path, c("NAME", est_vars), region = "county:*",
                      regionin = "state:*", cache = dat_cache)
    dat <- dat |>
      mutate(fips = fips5(state, county), acs_end = end) |>
      relocate(fips, acs_end)

    acs_raw[[ey]][[tb]]    <- dat
    acs_labels[[ey]][[tb]] <- labs
    resolved[[length(resolved) + 1L]] <- tibble(
      acs_end = end, table = tb, n_vars = length(est_vars),
      concept = labs$concept[1], n_counties = nrow(dat)
    )
    message(sprintf("    %-8s %3d vars  %4d counties", tb, length(est_vars), nrow(dat)))
  }
}

resolved_df <- dplyr::bind_rows(resolved)
readr::write_csv(resolved_df, file.path(dr_dir, "acs_variables_resolved.csv"))
message("\n  wrote data-raw/acs_variables_resolved.csv (", nrow(resolved_df), " table-vintages)")
invisible(acs_raw)
