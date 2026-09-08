# ---------------------------------------------------------------------------
# 00_config.R  --  Shared configuration for the county_turnout_panel build
# ---------------------------------------------------------------------------
# Sourced by every other script in data-raw/. Defines:
#   * the election-year -> data-vintage mapping (the single source of truth)
#   * cache paths and a small idempotent download helper
#   * thin wrappers around the Census Bureau API (needs CENSUS_API_KEY)
#
# Nothing here downloads anything or has side effects beyond creating the
# cache directory tree.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(jsonlite)
})

# --- paths -----------------------------------------------------------------

dr_dir     <- normalizePath(file.path("data-raw"), mustWork = TRUE)
cache_dir  <- file.path(dr_dir, "cache")
cache_sub  <- c("medsl", "cvap", "acs", "rucc", "geo", "decennial")
for (s in cache_sub) dir.create(file.path(cache_dir, s), recursive = TRUE, showWarnings = FALSE)

cache_path <- function(sub, ...) file.path(cache_dir, sub, ...)

# --- election-year <-> vintage map (SINGLE SOURCE OF TRUTH) ---------------
#
# acs_end      : end year of the ACS 5-year window (== the API vintage year)
# acs_start    : start year of that window (documentation only)
# cvap_label   : the "<start>-<end>" label used in the CVAP special-tab file names
# cvap_dir_year: sub-directory year on www2.census.gov/.../rdo/datasets/
# rucc_vintage : USDA ERS Rural-Urban Continuum Code release used for ruralurban_cc
# dec_vintage  : decennial census used for the fixed rural/urban population split
# dec_dataset  : Census API dataset id for that decennial urban/rural table
#
# Rationale for the ACS windows: the 5-year window whose END year is closest
# to the election (ties -> the window ending in the election year).  This is
# exactly the mapping given in the build request.

vintages <- tibble::tribble(
  ~election_year, ~acs_start, ~acs_end, ~cvap_label,  ~cvap_dir_year, ~rucc_vintage, ~dec_vintage, ~dec_dataset,
  2008L,          2006L,      2010L,    "2006-2010",   2010L,          2003L,         2010L,        "dec/sf1",
  2012L,          2008L,      2012L,    "2008-2012",   2012L,          2013L,         2010L,        "dec/sf1",
  2016L,          2012L,      2016L,    "2012-2016",   2016L,          2013L,         2010L,        "dec/sf1",
  2020L,          2016L,      2020L,    "2016-2020",   2020L,          2023L,         2020L,        "dec/dhc",
  2024L,          2020L,      2024L,    "2020-2024",   2024L,          2023L,         2020L,        "dec/dhc"
)

# RUCC vintage reuse (documented, not hidden):
#   * 2003 codes  -> 2008 election
#   * 2013 codes  -> 2012 and 2016 elections   (ERS did not revise between them)
#   * 2023 codes  -> 2020 and 2024 elections   (ERS did not revise between them)
# ERS revises the continuum roughly once per decade.

# --- generic cached download --------------------------------------------------

download_cached <- function(url, dest, mode = "wb", quiet = FALSE, force = FALSE,
                            headers = character()) {
  if (!force && file.exists(dest) && file.info(dest)$size > 0) {
    if (!quiet) message("  cached: ", basename(dest))
    return(invisible(dest))
  }
  if (!quiet) message("  downloading: ", url)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  req <- httr2::request(url) |>
    httr2::req_timeout(600) |>
    httr2::req_retry(max_tries = 4, backoff = ~ 2^.x) |>
    httr2::req_user_agent("pintervals data-raw build (https://github.com/Doktorandahl/pintervals)")
  if (length(headers)) req <- httr2::req_headers(req, !!!headers)
  resp <- httr2::req_perform(req, path = dest)
  httr2::resp_check_status(resp)
  invisible(dest)
}

# --- Census Bureau API ----------------------------------------------------------

census_key <- function() {
  k <- Sys.getenv("CENSUS_API_KEY")
  if (!nzchar(k))
    stop("CENSUS_API_KEY is not set. Add it to ~/.Renviron (see data-raw/README.md).",
         call. = FALSE)
  k
}

# One GET against the Census API, returned as a tibble (all columns character).
# `path` is e.g. "2016/acs/acs5" or "2010/dec/sf1"; `vars` a character vector of
# API variable names; `region` / `regionin` the geography predicates.
census_get <- function(path, vars, region = "county:*", regionin = NULL,
                       cache = NULL, force = FALSE) {
  if (!is.null(cache) && !force && file.exists(cache)) return(read_rds(cache))
  stopifnot(length(vars) >= 1)
  # API allows up to 50 variables per call; chunk and join on the geo columns.
  geo_cols <- if (grepl("county", region)) c("state", "county") else "state"
  chunks <- split(vars, ceiling(seq_along(vars) / 45L))
  out <- NULL
  for (ch in chunks) {
    q <- list(get = paste(c(ch), collapse = ","), key = census_key())
    q[["for"]] <- region
    if (!is.null(regionin)) q[["in"]] <- regionin
    url <- httr2::request(paste0("https://api.census.gov/data/", path)) |>
      httr2::req_url_query(!!!q) |>
      httr2::req_timeout(300) |>
      httr2::req_retry(max_tries = 4, backoff = ~ 2^.x)
    resp <- httr2::req_perform(url)
    httr2::resp_check_status(resp)
    m <- jsonlite::fromJSON(httr2::resp_body_string(resp))
    df <- as_tibble(as.data.frame(m[-1, , drop = FALSE], stringsAsFactors = FALSE))
    names(df) <- m[1, ]
    out <- if (is.null(out)) df else left_join(out, df, by = intersect(names(out), names(df)))
  }
  if (!is.null(cache)) write_rds(out, cache, compress = "xz")
  out
}

# Variable-level metadata for one API table group (used to pick cells by label,
# so the transforms stay robust to table-structure differences across vintages).
census_group_vars <- function(path, group, cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) return(read_rds(cache))
  url <- sprintf("https://api.census.gov/data/%s/groups/%s.json", path, group)
  j <- jsonlite::fromJSON(url)
  v <- j$variables
  tb <- tibble(
    name  = names(v),
    label = map_chr(v, "label", .default = NA_character_),
    concept = map_chr(v, "concept", .default = NA_character_)
  ) |>
    filter(str_ends(name, "E"), !name %in% c("NAME")) |>
    arrange(name)
  if (!is.null(cache)) write_rds(tb, cache, compress = "xz")
  tb
}

fips5 <- function(state, county) paste0(str_pad(state, 2, "left", "0"),
                                        str_pad(county, 3, "left", "0"))

message("config loaded: ", nrow(vintages), " election waves ",
        paste(vintages$election_year, collapse = ", "))
