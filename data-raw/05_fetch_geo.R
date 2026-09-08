# ---------------------------------------------------------------------------
# 05_fetch_geo.R -- County centroids (longitude / latitude)
# ---------------------------------------------------------------------------
# Source : Census Bureau Gazetteer Files, county-level "internal point"
#          coordinates (INTPTLAT / INTPTLONG). U.S. Government work, public
#          domain.
# Output : data-raw/cache/geo/<year>_Gaz_counties_national.zip (+ extracted .txt)
#
# Gazetteer vintage per election wave. Gazetteer files begin in 2012, so the
# 2008 wave borrows the 2012 file (county internal points move < 1 km between
# vintages when boundaries are unchanged; documented in the codebook).
#   2008 -> 2012 gaz | 2012 -> 2013 gaz | 2016 -> 2016 gaz
#   2020 -> 2020 gaz | 2024 -> 2024 gaz
#
# region / division / geo_group are NOT downloaded here: they are a fixed
# state-level classification inherited from the existing bundled `county_turnout`
# dataset (U.S. Census Bureau regions & divisions; `geo_group` is a
# pintervals-specific coarse grouping). See 20_assemble_panel.R.
# ---------------------------------------------------------------------------

if (!exists("vintages")) source(file.path("data-raw", "00_config.R"))

gaz_year <- c(`2008` = 2012L, `2012` = 2013L, `2016` = 2016L,
              `2020` = 2020L, `2024` = 2024L)

gaz_base <- "https://www2.census.gov/geo/docs/maps-data/data/gazetteer"

geo_files <- character()
for (ey in names(gaz_year)) {
  gy   <- gaz_year[[ey]]
  zip  <- cache_path("geo", sprintf("%d_Gaz_counties_national.zip", gy))
  txt  <- cache_path("geo", sprintf("%d_Gaz_counties_national.txt", gy))
  download_cached(sprintf("%s/%d_Gazetteer/%d_Gaz_counties_national.zip", gaz_base, gy, gy), zip)
  if (!file.exists(txt)) {
    nm <- utils::unzip(zip, list = TRUE)$Name
    utils::unzip(zip, files = nm[grepl("counties", nm)], exdir = cache_path("geo"), junkpaths = TRUE)
    # normalise extracted name
    got <- list.files(cache_path("geo"), pattern = sprintf("(?i)%d.*counties.*\\.txt$", gy), full.names = TRUE)
    if (length(got) && !file.exists(txt)) file.rename(got[1], txt)
  }
  geo_files[ey] <- txt
  message(sprintf("  gazetteer %s wave <- %d file", ey, gy))
}

.chk <- readr::read_tsv(geo_files[["2024"]], show_col_types = FALSE, n_max = 3,
                        name_repair = "minimal")
message("  gazetteer columns: ", paste(trimws(names(.chk)), collapse = ", "))
invisible(geo_files)
