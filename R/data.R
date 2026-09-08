#' Election-Year Democracy Indicators from V-Dem (1946–2024)
#'
#' A sample of election years from the V-Dem dataset covering 2,680 country-years between 1946 and 2024. Includes a range of democracy indices and related variables measured during years in which national elections were held.
#'
#' @format ## ´elections´ A tibble with 2,680 rows and 21 variables:
#' \describe{
#'   \item{country_name}{Country name}
#'   \item{year}{Election year}
#'   \item{v2x_polyarchy}{Electoral democracy index}
#'   \item{v2x_libdem}{Liberal democracy index}
#'   \item{v2x_partipdem}{Participatory democracy index}
#'   \item{v2x_delibdem}{Deliberative democracy index}
#'   \item{v2x_egaldem}{Egalitarian democracy index}
#'   \item{v2xel_frefair}{Free and fair elections index}
#'   \item{v2x_frassoc_thick}{Freedom of association index}
#'   \item{v2x_elecoff}{Elected officials index}
#'   \item{v2eltrnout}{Voter turnout (V-Dem)}
#'   \item{v2x_accountability}{Accountability index}
#'   \item{v2xps_party}{Party system institutionalization}
#'   \item{v2x_civlib}{Civil liberties index}
#'   \item{v2x_corr}{Control of corruption index}
#'   \item{v2x_rule}{Rule of law index}
#'   \item{v2x_neopat}{Neo-patrimonial rule index}
#'   \item{v2x_suffr}{Suffrage index}
#'   \item{turnout}{Turnout percentage (external source)}
#'   \item{hog_lost}{Factor indicating if head of government lost election}
#'   \item{hog_lost_num}{Numeric version of \code{hog_lost}}
#' }
#'
#' @source Data derived from the Varieties of Democracy (V-Dem) dataset, version 15, filtered to election years between 1946 and 2024. <https://www.v-dem.net/data/the-v-dem-dataset/>
"elections"

#' U.S. County-Level Turnout and Demographic Context (MIT Election Lab 2018 Election Analysis Dataset + Additions)
#'
#' @description
#' A county-level dataset (U.S.) with voter turnout and sociodemographic covariates.
#'
#' @format A tibble with 3,107 rows and 22 variables:
#' \describe{
#'   \item{state}{State name.}
#'   \item{county}{County name.}
#'   \item{fips}{County FIPS code.}
#'   \item{turnout}{Observed turnout (proportion). Total presidential votes cast divided by the county Citizen Voting-Age Population (CVAP) from the 2012--2016 ACS special tabulation. (Earlier versions of this help page described the denominator as total population; that was incorrect -- see NEWS for 1.3.0.)}
#'   \item{total_population}{Total county population.}
#'   \item{nonwhite_pct}{Percent non-white population.}
#'   \item{foreignborn_pct}{Percent foreign-born population.}
#'   \item{female_pct}{Percent female population.}
#'   \item{age29andunder_pct}{Percent of population aged 29 or under.}
#'   \item{age65andolder_pct}{Percent of population aged 65 or older.}
#'   \item{median_hh_inc}{Median household income.}
#'   \item{clf_unemploy_pct}{Percent unemployed in the civilian labor force.}
#'   \item{lesscollege_pct}{Percent with less than college education.}
#'   \item{lesshs_pct}{Percent with less than high school education.}
#'   \item{rural_pct}{Percent rural.}
#'   \item{ruralurban_cc}{Rural–urban continuum code.}
#'   \item{predicted_turnout}{LOO-CV random-forest prediction of `turnout` (see Details).}
#'   \item{division}{U.S. Census division.}
#'   \item{region}{U.S. Census region.}
#'   \item{geo_group}{Additional coarse geographic grouping variable (added).}
#'   \item{longitude}{County centroid longitude (added).}
#'   \item{latitude}{County centroid latitude (added).}
#' }
#'
#' @details
#' \strong{Superseded.} This single-election (2016) cross-section is retained for
#' backward compatibility and is still used by the package vignette and examples.
#' For new work spanning multiple elections, prefer \code{\link{county_turnout_panel}},
#' a five-wave (2008--2024) county panel whose covariates are pulled directly from
#' the U.S. Census Bureau API. \code{county_turnout_panel$turnout} uses the same
#' votes / Citizen Voting-Age Population definition as \code{county_turnout$turnout}
#' and reproduces its 2016 values for all but a handful of counties (see NEWS 1.3.0).
#'
#' The dataset is based on the MIT Election Lab "2018 Election Analysis dataset" file, with four additions: (1) `turnout`, calculated as the number of presidential votes cast divided by the Citizen Voting-Age Population, (2) `geo_group`, a coarse geographic grouping variable for the counties, (3) county centroid coordinates (`longitude`, `latitude`), and (4) `predicted_turnout`.
#' The variable `predicted_turnout` is generated using leave-one-out
#' cross-validation (LOO-CV). For each county a random forest
#' model is fit on the remaining counties with `turnout` as the outcome
#' and all available *non-geographic* covariates as predictors. The fitted model
#' is then used to predict turnout for the held-out county. Geographic
#' features are excluded from the predictor set to avoid leaking spatial
#' information into the prediction target. Concretely, identifiers and geographic variables (e.g., `state`, `county`, `fips`, `division`, `region`, `geo_group`, `longitude`, `latitude`) are excluded from the predictor set.
#'
#' Below is example code (using `foreach`) to reproduce `predicted_turnout`. This is computationally expensive for LOO-CV; parallel execution is recommended.
#'
#' library(dplyr)
#' library(ranger)
#' library(foreach)
#' library(pintervals)
#'
#' data("county_turnout", package = "pintervals")
#'
#' # Choose predictors: all numeric covariates except turnout + geographic/id vars
#' model_data <- county_turnout \code{\%>\%}
#' select(-c(state, county, fips, division, region, geo_group, longitude, latitude))
#'
#' set.seed(101010) # The meaning of life in binary
#'
#' pred_loo <- foreach(.i = seq_len(nrow(county_turnout)),
#'                     .final = unlist) \%do\% \{
#'
#' train <- model_data[-.i, , drop = FALSE]
#' test  <- model_data[ .i, , drop = FALSE]
#'
#'   fit <- ranger(
#'     formula = turnout ~ .,
#'     data = train
#'   )
#'
#'   predict(fit, data = test)$predictions[[1]]
#'
#' \}
#'
#'
#' county_turnout <- county_turnout \code{\%>\%} mutate(predicted_turnout = pred_loo)
#'
#'
#' @source
#' The base covariates originate from the MEDSL "2018 election context" file:
#' \url{https://github.com/MEDSL/2018-elections-unoffical/blob/master/election-context-2018.md}.
#' The variables `geo_group`, `longitude`, `latitude`, and `predicted_turnout`
#' are additions.
#'
#' @usage data(county_turnout)
#'
#' @keywords datasets
"county_turnout"

#' U.S. County-Level Presidential Turnout and Demographic Context, 2008-2024 Panel
#'
#' @description
#' A tidy panel of U.S. county-level presidential-election voter turnout and
#' sociodemographic context for the five presidential elections of 2008, 2012,
#' 2016, 2020 and 2024. One row per county (county-equivalent) per election year.
#' This is the multi-election generalisation of \code{\link{county_turnout}}; every covariate
#' is re-derived directly from the U.S. Census Bureau API.
#'
#' @format A tibble with 15,543 rows and 22 variables:
#' \describe{
#'   \item{state}{State name (or "District of Columbia").}
#'   \item{county}{County / county-equivalent name (trailing " County"/" Parish" removed). Not unique on its own; pair with \code{fips}.}
#'   \item{fips}{5-character county FIPS code, stored as a character string with leading zeros (e.g. \code{"01001"}). (In \code{county_turnout} this column is numeric.)}
#'   \item{election_year}{Presidential election year: 2008, 2012, 2016, 2020 or 2024. Panel key together with \code{fips}.}
#'   \item{turnout}{Turnout rate = total presidential ballots cast (MEDSL) / Citizen Voting-Age Population (CVAP special tabulation), using the same definition as \code{county_turnout$turnout}. Always \eqn{\le 1} where present: 25 micro-county-years whose recorded ballots exceeded the noisy 5-year ACS CVAP estimate (turnout > 1) have been \strong{dropped} from the panel. \code{NA} for 14 rows (Kalawao Co., HI in every wave; the 9 Connecticut planning regions in 2024) - these are \strong{retained}.}
#'   \item{total_population}{Total resident population (ACS table B01003), integer.}
#'   \item{nonwhite_pct}{Percent not non-Hispanic White (from ACS B03002).}
#'   \item{foreignborn_pct}{Percent foreign-born (ACS B05002).}
#'   \item{female_pct}{Percent female (ACS B01001).}
#'   \item{age29andunder_pct}{Percent aged 29 or younger (ACS B01001).}
#'   \item{age65andolder_pct}{Percent aged 65 or older (ACS B01001).}
#'   \item{median_hh_inc}{Median household income, inflation-adjusted to the ACS window's final year (ACS B19013), integer. \code{NA} for 2 Census-suppressed micro-county rows.}
#'   \item{clf_unemploy_pct}{Unemployed as a percent of the civilian labor force (ACS B23025; B23001 for the 2008 wave).}
#'   \item{lesscollege_pct}{Percent of the population 25+ with less than a bachelor's degree (ACS B15003; B15002 for the 2008 wave).}
#'   \item{lesshs_pct}{Percent of the population 25+ with less than a high-school diploma (ACS B15003; B15002 for the 2008 wave).}
#'   \item{rural_pct}{Percent of population classified rural by the nearest decennial census (table P2): 2010 census for the 2008/2012/2016 waves, 2020 census for the 2020/2024 waves. Held fixed across waves that share a census. \code{NA} for the 9 Connecticut planning regions in 2024.}
#'   \item{ruralurban_cc}{USDA ERS Rural-Urban Continuum Code (integer 1-9; lower = more metropolitan). Vintage reused across adjacent elections: 2003 -> 2008; 2013 -> 2012 & 2016; 2023 -> 2020 & 2024.}
#'   \item{division}{U.S. Census Bureau division (factor, 9 levels). State-level; inherited from \code{county_turnout}.}
#'   \item{region}{U.S. Census Bureau region (factor, 4 levels). State-level; inherited from \code{county_turnout}.}
#'   \item{geo_group}{Coarse \code{pintervals}-specific state grouping (factor, 15 levels). Inherited from \code{county_turnout}.}
#'   \item{longitude}{County internal-point (approximate centroid) longitude, decimal degrees (Census Gazetteer).}
#'   \item{latitude}{County internal-point latitude, decimal degrees (Census Gazetteer).}
#' }
#'
#' @details
#' Each covariate is matched to the ACS 5-year vintage whose window ends nearest
#' the election: 2006-2010 (2008), 2008-2012 (2012), 2012-2016 (2016),
#' 2016-2020 (2020), 2020-2024 (2024).
#'
#' \strong{Coverage.} 50 states + DC. \strong{Alaska is excluded} - MEDSL reports
#' Alaska by state-house district, not by borough, so no borough-level
#' presidential returns exist (the original \code{county_turnout} also omits
#' Alaska). Island territories are excluded (no presidential electors).
#'
#' \strong{ACS covariates are pulled directly from the U.S. Census Bureau API,
#' not via IPUMS NHGIS.} NHGIS's terms of use prohibit redistributing extracted
#' data outside narrow journal exceptions, so ACS data obtained through NHGIS
#' cannot be bundled in a CRAN package; the identical ACS detailed tables from
#' the Census Bureau's own API are a U.S. Government work in the public domain.
#' The reconstruction was validated against \code{county_turnout} (which is built
#' from an NHGIS-sourced file): for the 2012-2016 ACS vintage all ten shared
#' covariates match to floating-point tolerance across all 3,107 shared counties.
#'
#' \strong{Geography / FIPS harmonisation.} Geography is kept native to each wave
#' (the matched ACS county list is the spine). Documented handling of known
#' changes in the window:
#' \itemize{
#'   \item \emph{Shannon County, SD (46113) -> Oglala Lakota County, SD (46102), 2015}: sources still using 46113 (MEDSL, RUCC 2013, 2010 census) are recoded to 46102 for the 2016/2020/2024 waves; 2008/2012 use 46113.
#'   \item \emph{Bedford (independent) city, VA (51515) -> Bedford County, VA (51019), 2013}: separate row in 2008/2012 (native); folded into 51019 for 2016; 51019 only thereafter.
#'   \item \emph{Valdez-Cordova Census Area, AK split (2019)} and \emph{Wade Hampton -> Kusilvak, AK (2015)}: Alaska is excluded, so these are noted only for completeness.
#'   \item \emph{Connecticut: 8 counties -> 9 planning regions (2022-vintage products)}: kept native per wave - the 2008-2020 waves use the 8 legacy counties (09001-09015), the 2024 wave uses the 9 planning regions (09110-09190). Connecticut is therefore not a consistent 5-wave panel: 2024 \code{turnout} and \code{rural_pct} are \code{NA} for the planning regions, and CT's 2020 \code{ruralurban_cc} is taken from RUCC 2013 (legacy counties).
#'   \item \emph{Kansas City, MO} reports as its own unit (FIPS 2938000; mis-coded 36000 in the 2024 MEDSL file); its ballots are redistributed to Cass/Clay/Jackson/Platte counties by each county's share of the Kansas City place population from the nearest decennial census.
#' }
#' Native geography per wave is 3114 / 3114 / 3113 / 3113 / 3114 county-equivalents
#' for 2008 / 2012 / 2016 / 2020 / 2024 (the Bedford (2016) and Connecticut (2024)
#' changes above); after dropping the 6 / 3 / 4 / 6 / 6 county-years with
#' turnout > 1, the shipped per-year row counts are
#' 3108 / 3111 / 3109 / 3107 / 3108.
#'
#' \strong{Rurality-code vintage caveat.} USDA ERS revises the Rural-Urban
#' Continuum roughly once per decade, so \code{ruralurban_cc} is shared by the
#' 2012 & 2016 waves (2013 codes) and by the 2020 & 2024 waves (2023 codes).
#'
#' A full variable-by-variable codebook is installed with the package at
#' \code{system.file("extdata", "county_turnout_panel_codebook.md", package = "pintervals")}.
#' The dataset is rebuilt end-to-end by \code{data-raw/county_turnout_panel.R};
#' see \code{data-raw/README.md} for the exact upstream citations and the
#' regeneration procedure.
#'
#' @source
#' \itemize{
#'   \item \strong{Turnout numerator}: MEDSL, "County Presidential Election Returns 2000-2024", Harvard Dataverse, \doi{10.7910/DVN/VOQCHQ}. License: \strong{CC0 1.0}.
#'   \item \strong{Turnout denominator}: U.S. Census Bureau, Citizen Voting-Age Population (CVAP) Special Tabulation from the ACS, \url{https://www.census.gov/programs-surveys/decennial-census/about/voting-rights/cvap.html}. \strong{Public domain} (U.S. Government work).
#'   \item \strong{Demographic covariates}: U.S. Census Bureau, American Community Survey 5-year estimates (detailed tables B01001, B01003, B03002, B05002, B15002/B15003, B19013, B23001/B23025), obtained directly from the Census Bureau API (\url{https://www.census.gov/data/developers/data-sets/acs-5year.html}). \strong{Public domain}. Deliberately not obtained via IPUMS NHGIS, to comply with NHGIS's non-redistribution terms.
#'   \item \strong{rural_pct}: U.S. Census Bureau, decennial census table P2 (Urban and Rural), 2010 and 2020, via the Census Bureau API. \strong{Public domain}.
#'   \item \strong{ruralurban_cc}: USDA Economic Research Service, Rural-Urban Continuum Codes (2003, 2013, 2023 releases), \url{https://www.ers.usda.gov/data-products/rural-urban-continuum-codes/}. \strong{Public domain} (U.S. Government work).
#'   \item \strong{longitude / latitude}: U.S. Census Bureau Gazetteer Files, \url{https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.html}. \strong{Public domain}.
#'   \item \strong{division / region / geo_group}: inherited verbatim from \code{pintervals::county_turnout} (U.S. Census Bureau regions and divisions; \code{geo_group} is a \code{pintervals} addition).
#' }
#'
#' @seealso \code{\link{county_turnout}} for the original single-election cross-section.
#'
#' @usage data(county_turnout_panel)
#'
#' @keywords datasets
"county_turnout_panel"
