#' Election-year democracy indicators from V-Dem (1946–2024)
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
