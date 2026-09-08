# ---------------------------------------------------------------------------
# 04_fetch_rucc.R -- USDA ERS Rural-Urban Continuum Codes
# ---------------------------------------------------------------------------
# Source : USDA Economic Research Service, Rural-Urban Continuum Codes.
#          U.S. Government work, public domain.
# Access : direct file download from ers.usda.gov.
# Output : data-raw/cache/rucc/rucc_<vintage>.<ext>
#
# Vintages used (ERS revises ~once per decade; reuse is documented, not hidden):
#   2003 codes -> 2008 election
#   2013 codes -> 2012 & 2016 elections
#   2023 codes -> 2020 & 2024 elections
# ---------------------------------------------------------------------------

if (!exists("vintages")) source(file.path("data-raw", "00_config.R"))

# Stable ERS media endpoints (the "?v=" cache-buster is required by the CDN).
rucc_sources <- tibble::tribble(
  ~vintage, ~url,                                                                            ~ext,
  2003L,    "https://www.ers.usda.gov/media/5770/2003-rural-urban-continuum-codes.xls?v=94819",  "xls",
  2013L,    "https://www.ers.usda.gov/media/5769/2013-rural-urban-continuum-codes.xls?v=88984",  "xls",
  2023L,    "https://www.ers.usda.gov/media/5768/2023-rural-urban-continuum-codes.csv?v=87920",  "csv"
)

rucc_files <- character()
for (i in seq_len(nrow(rucc_sources))) {
  s    <- rucc_sources[i, ]
  dest <- cache_path("rucc", sprintf("rucc_%d.%s", s$vintage, s$ext))
  download_cached(s$url, dest)
  rucc_files[as.character(s$vintage)] <- dest
}

message("  RUCC vintages cached: ", paste(names(rucc_files), collapse = ", "))
invisible(rucc_files)
