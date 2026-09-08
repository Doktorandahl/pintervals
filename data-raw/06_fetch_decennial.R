# ---------------------------------------------------------------------------
# 06_fetch_decennial.R -- Rural / urban population split (for rural_pct)
# ---------------------------------------------------------------------------
# The ACS contains no rural/urban classification, so rural_pct is taken from
# the nearest decennial census "Urban and Rural" table (P2) and held fixed
# across the waves that share a census:
#   2008 / 2012 / 2016 waves  <- 2010 Census (dec/sf1  table P2:  P002001/P002005)
#   2020 / 2024 waves         <- 2020 Census (dec/dhc  table P2:  P2_001N/P2_003N)
#
# Source : Census Bureau decennial census, U.S. Government work, public domain.
# Access : Census Bureau API (needs CENSUS_API_KEY).
# Output : data-raw/cache/decennial/rural_<censusYear>.rds
# ---------------------------------------------------------------------------

if (!exists("vintages")) source(file.path("data-raw", "00_config.R"))

fetch_rural <- function(census_year, dataset) {
  dest <- cache_path("decennial", sprintf("rural_%d.rds", census_year))
  if (file.exists(dest)) { message("  cached: ", basename(dest)); return(read_rds(dest)) }

  if (dataset == "dec/sf1") {
    vars <- c(total = "P002001", urban = "P002002", rural = "P002005")
  } else {                       # dec/dhc (2020)
    vars <- c(total = "P2_001N", urban = "P2_002N", rural = "P2_003N")
  }
  raw <- census_get(sprintf("%d/%s", census_year, dataset), unname(vars),
                    region = "county:*", regionin = "state:*")
  out <- raw |>
    rename(!!!vars) |>                       # vars = c(total="P002001", ...)
    mutate(across(all_of(names(vars)), as.numeric),
           fips = fips5(state, county),
           rural_pct = 100 * rural / total,
           census_year = census_year) |>
    select(fips, census_year, pop_total_dec = total, pop_rural_dec = rural, rural_pct)
  write_rds(out, dest, compress = "xz")
  message(sprintf("  decennial %d rural: %d counties", census_year, nrow(out)))
  out
}

rural_2010 <- fetch_rural(2010L, "dec/sf1")
rural_2020 <- fetch_rural(2020L, "dec/dhc")
rural_by_census <- list(`2010` = rural_2010, `2020` = rural_2020)
invisible(rural_by_census)
