# ---------------------------------------------------------------------------
# 07_fetch_kc_split.R -- Kansas City, MO -> constituent-county weights
# ---------------------------------------------------------------------------
# MEDSL reports Kansas City, MO as its own unit (FIPS 2938000 in 2000-2020;
# mis-coded as 3600000 / "KANSAS CITY" in the 2024 file -- an upstream data
# error, flagged in the build notes). Kansas City straddles four counties
# (Cass, Clay, Jackson, Platte) whose own MEDSL rows EXCLUDE the Kansas City
# vote. To get a usable county turnout numerator, the Kansas City ballots are
# redistributed to those four counties in proportion to the share of the
# Kansas City *place* population that each county contains, from the nearest
# decennial census (2010 -> 2008/2012/2016 waves; 2020 -> 2020/2024 waves).
#
# Source : Census Bureau decennial census, "place (or part)" within county.
# Output : data-raw/cache/decennial/kc_split_<censusYear>.rds
# ---------------------------------------------------------------------------

if (!exists("vintages")) source(file.path("data-raw", "00_config.R"))

fetch_kc_split <- function(census_year) {
  dest <- cache_path("decennial", sprintf("kc_split_%d.rds", census_year))
  if (file.exists(dest)) { message("  cached: ", basename(dest)); return(read_rds(dest)) }

  if (census_year == 2010L) { ds <- "2010/dec/sf1"; popvar <- "P001001" }
  else                      { ds <- "2020/dec/dhc"; popvar <- "P1_001N" }

  resp <- httr2::request("https://api.census.gov/data") |>
    httr2::req_url_path_append(strsplit(ds, "/")[[1]]) |>
    httr2::req_url_query(get = paste0("NAME,", popvar),
                         `for` = "county (or part):*",
                         `in`  = "state:29 place:38000",
                         key   = census_key()) |>
    httr2::req_perform()
  m <- jsonlite::fromJSON(httr2::resp_body_string(resp))
  df <- as_tibble(as.data.frame(m[-1, , drop = FALSE], stringsAsFactors = FALSE))
  names(df) <- m[1, ]

  out <- df |>
    transmute(fips = fips5("29", .data[["county (or part)"]]),
              kc_pop = as.numeric(.data[[popvar]])) |>
    mutate(kc_weight = kc_pop / sum(kc_pop), census_year = census_year)
  write_rds(out, dest, compress = "xz")
  message(sprintf("  KC split %d: %s", census_year,
                  paste(sprintf("%s=%.3f", out$fips, out$kc_weight), collapse = " ")))
  out
}

kc_split <- list(`2010` = fetch_kc_split(2010L), `2020` = fetch_kc_split(2020L))
invisible(kc_split)
