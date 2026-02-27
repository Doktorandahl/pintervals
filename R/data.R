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
#'   \item{turnout}{Observed turnout (proportion). Calculated as total votes cast divided by total population (not voting-age population).}
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
#' The dataset is based on the MIT Election Lab "2018 Election Analysis dataset" file, with four additions: (1) `turnout`, calculated as the number of votes cast divided by the total population, (2) `geo_group`, a coarse geographic grouping variable for the counties, (3) county centroid coordinates (`longitude`, `latitude`), and (4) `predicted_turnout`.
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
#' dat <- county_turnout  # replace with your object name
#'
#' # Choose predictors: all numeric covariates except turnout + geographic/id vars
#' dat2 <- dat \code{|>}
#'  select(-c(state, county, fips, division, region, geo_group, longitude, latitude))
#'
#' set.seed(101010) # The meaning of life in binary
#'
#' pred_loo <- foreach(.i = seq_len(nrow(dat)),
#'                     .final = unlist) %do% \{
#'
#' train <- dat2[-.i, , drop = FALSE]
#'   test  <- dat2[ .i, , drop = FALSE]
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
#' dat <- dat \code{|>} mutate(predicted_turnout = pred_loo)
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
