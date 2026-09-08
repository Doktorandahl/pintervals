# ---------------------------------------------------------------------------
# 01_fetch_medsl.R -- Vote totals (turnout numerator)
# ---------------------------------------------------------------------------
# Source : MEDSL "County Presidential Election Returns 2000-2024"
#          Harvard Dataverse, doi:10.7910/DVN/VOQCHQ   (License: CC0 1.0)
# Access : Dataverse data-access API (NOT the web UI).
#
# NOTE ON ACCESS: the dataset content is CC0, but Harvard Dataverse attaches a
# *required* download guestbook (name / email / institution / position) to it.
# This is an access formality, not a licence term. It is satisfied here by
# POSTing a guestbook response to the access endpoint, which returns a short-
# lived signed S3 URL. No Dataverse account or API token is required. Edit the
# `gb_response` list below if you would rather record your own details.
#
# Output : data-raw/cache/medsl/countypres_2000-2024.tab
#          data-raw/cache/medsl/sources-president.tab      (per-state provenance)
#          data-raw/cache/medsl/dataset-metadata.json      (version + license)
# ---------------------------------------------------------------------------

if (!exists("vintages")) source(file.path("data-raw", "00_config.R"))

medsl_doi <- "doi:10.7910/DVN/VOQCHQ"
dv_base   <- "https://dataverse.harvard.edu/api"

gb_response <- list(
  name        = "pintervals package data-raw build (automated)",
  email       = Sys.getenv("DATAVERSE_GUESTBOOK_EMAIL", "david@randahl.org"),
  institution = "Uppsala University, Department of Peace and Conflict Research",
  position    = "Researcher"
)

## 1. dataset metadata (pin the version, record the license) -----------------
meta_dest <- cache_path("medsl", "dataset-metadata.json")
download_cached(
  sprintf("%s/datasets/:persistentId/?persistentId=%s", dv_base, medsl_doi),
  meta_dest, force = TRUE
)
meta  <- jsonlite::fromJSON(meta_dest)
lv    <- meta$data$latestVersion
files <- lv$files
message(sprintf("  MEDSL dataset version %s.%s, license: %s",
                lv$versionNumber, lv$versionMinorNumber, lv$license$name))
stopifnot(identical(toupper(lv$license$rightsIdentifier), "CC0-1.0"))

## 2. resolve file ids by label (ids are version-specific) ------------------
want <- c(returns = "countypres_2000-2024.tab",
          sources = "sources-president.tab")
ids  <- setNames(files$dataFile$id[match(want, files$label)], names(want))
stopifnot(!anyNA(ids))

## 3. download each file through a guestbook-response signed URL ------------
medsl_signed_url <- function(file_id) {
  resp <- httr2::request(dv_base) |>
    httr2::req_url_path_append("access", "datafile", file_id) |>
    httr2::req_body_json(list(guestbookResponse = gb_response)) |>
    httr2::req_retry(max_tries = 3) |>
    httr2::req_perform()
  out <- jsonlite::fromJSON(httr2::resp_body_string(resp))
  if (!identical(out$status, "OK") || is.null(out$data$signedUrl))
    stop("Dataverse did not return a signed URL for file ", file_id, call. = FALSE)
  out$data$signedUrl
}

for (nm in names(ids)) {
  dest <- cache_path("medsl", want[[nm]])
  if (file.exists(dest) && file.info(dest)$size > 0) {
    message("  cached: ", basename(dest)); next
  }
  message("  downloading (guestbook): ", want[[nm]])
  download_cached(medsl_signed_url(ids[[nm]]), dest)
}

## 4. sanity read ---------------------------------------------------------------
returns_raw <- readr::read_tsv(cache_path("medsl", want[["returns"]]),
                               show_col_types = FALSE, guess_max = 1e6)
message(sprintf("  MEDSL returns: %s rows, years %s",
                format(nrow(returns_raw), big.mark = ","),
                paste(range(returns_raw$year), collapse = "-")))
stopifnot(all(c(2008, 2012, 2016, 2020, 2024) %in% returns_raw$year))
attr(returns_raw, "medsl_version") <- sprintf("%s.%s", lv$versionNumber, lv$versionMinorNumber)
invisible(returns_raw)
