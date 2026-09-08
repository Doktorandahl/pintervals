# ---------------------------------------------------------------------------
# 02_fetch_cvap.R -- Citizen Voting-Age Population (turnout denominator)
# ---------------------------------------------------------------------------
# Source : Census Bureau CVAP Special Tabulation from the ACS (via the
#          Redistricting Data Office, census.gov/rdo). U.S. Government work,
#          public domain.
# Access : direct file download of the published CSV zip per ACS vintage.
# Output : data-raw/cache/cvap/CVAP_<label>_ACS_csv_files.zip  (+ extracted)
#
# One vintage per election wave (see 00_config.R $vintages):
#   2006-2010 -> 2008 | 2008-2012 -> 2012 | 2012-2016 -> 2016
#   2016-2020 -> 2020 | 2020-2024 -> 2024
# ---------------------------------------------------------------------------

if (!exists("vintages")) source(file.path("data-raw", "00_config.R"))

rdo_base <- "https://www2.census.gov/programs-surveys/decennial/rdo/datasets"

cvap_url <- function(dir_year, label) {
  # www2.census.gov/.../rdo/datasets/<dir_year>/<dir_year>-cvap/CVAP_<label>_ACS_csv_files.zip
  sprintf("%s/%d/%d-cvap/CVAP_%s_ACS_csv_files.zip", rdo_base, dir_year, dir_year, label)
}

cvap_files <- character()

for (i in seq_len(nrow(vintages))) {
  v      <- vintages[i, ]
  label  <- v$cvap_label
  zip    <- cache_path("cvap", sprintf("CVAP_%s_ACS_csv_files.zip", label))
  exdir  <- cache_path("cvap", label)

  download_cached(cvap_url(v$cvap_dir_year, label), zip)

  csv <- list.files(exdir, pattern = "(?i)county\\.csv$", full.names = TRUE, recursive = TRUE)
  if (length(csv) == 0) {
    dir.create(exdir, showWarnings = FALSE)
    # extract ONLY the county-level CSV -- the block-group / tract members of
    # these zips are hundreds of MB each and are not used here.
    members <- utils::unzip(zip, list = TRUE)$Name
    want_member <- members[grepl("(?i)county\\.csv$", members)]
    if (length(want_member) == 0)
      stop("No county-level CVAP CSV inside ", basename(zip), call. = FALSE)
    utils::unzip(zip, files = want_member[1], exdir = exdir, junkpaths = TRUE)
    csv <- list.files(exdir, pattern = "(?i)county\\.csv$", full.names = TRUE, recursive = TRUE)
  }
  cvap_files[label] <- csv[1]
  message(sprintf("  CVAP %s: %s", label, basename(csv[1])))
}

# quick structural check on the most recent vintage
.chk <- readr::read_csv(cvap_files[[nrow(vintages)]], show_col_types = FALSE, n_max = 5)
message("  CVAP county CSV columns: ", paste(names(.chk), collapse = ", "))
invisible(cvap_files)
