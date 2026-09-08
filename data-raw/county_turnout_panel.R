# ===========================================================================
# county_turnout_panel.R  --  end-to-end build of data/county_turnout_panel.rda
# ===========================================================================
# Idempotent: raw downloads are cached under data-raw/cache/ (git-ignored) and
# re-used on subsequent runs. Delete a cache sub-directory to force a refresh.
#
# Requires:  CENSUS_API_KEY in ~/.Renviron  (Census Bureau API)
# Run from the package root:  Rscript data-raw/county_turnout_panel.R
#
# Stages (see the individual scripts for source URLs and licences):
#   00  config: election-year -> data-vintage map, helpers
#   01  fetch MEDSL county presidential returns   (Harvard Dataverse, CC0)
#   02  fetch CVAP special tabulations            (census.gov/rdo, public domain)
#   03  fetch ACS 5-year detailed tables          (Census API, public domain)
#   04  fetch USDA ERS Rural-Urban Continuum Codes (ers.usda.gov, public domain)
#   05  fetch Census Gazetteer county centroids   (census.gov, public domain)
#   06  fetch decennial urban/rural split         (Census API, public domain)
#   07  fetch Kansas City place-by-county weights  (Census API, public domain)
#   10  transform each source into a tidy county-year frame
#   20  assemble the panel (joins + FIPS harmonisation)
#   30  validate (schema, duplicates, turnout range, county counts)
#   -- write data/county_turnout_panel.rda + refresh the codebook
# ===========================================================================

t0 <- Sys.time()
stopifnot(file.exists("DESCRIPTION"))          # must run from package root

source("data-raw/00_config.R")
source("data-raw/01_fetch_medsl.R")
source("data-raw/02_fetch_cvap.R")
source("data-raw/03_fetch_acs.R")
source("data-raw/04_fetch_rucc.R")
source("data-raw/05_fetch_geo.R")
source("data-raw/06_fetch_decennial.R")
source("data-raw/07_fetch_kc_split.R")
source("data-raw/10_transform.R")
source("data-raw/20_assemble_panel.R")
source("data-raw/30_validate.R")

## --- write the dataset ----------------------------------------------------
county_turnout_panel <- tibble::as_tibble(county_turnout_panel)
save(county_turnout_panel, file = "data/county_turnout_panel.rda",
     compress = "xz", version = 3)
message(sprintf("wrote data/county_turnout_panel.rda (%.0f KB)",
                file.size("data/county_turnout_panel.rda") / 1024))

## The variable-by-variable codebook is hand-maintained at
## inst/extdata/county_turnout_panel_codebook.md -- keep it in sync when the
## schema or a source changes.

message(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
message("reminder: keep inst/extdata/county_turnout_panel_codebook.md in sync")
