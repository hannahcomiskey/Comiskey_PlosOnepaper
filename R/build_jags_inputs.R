#' Build JAGS metadata objects
#'
#' @description
#' From an indexed FP source data frame, build common metadata objects required
#' for JAGS / modeling: lookups, counts of subnational units per country,
#' match vectors used when indexing into parameter arrays, sequences of years, etc.
#'
#' @param df Data frame containing index columns produced by `add_index_variables()`.
#' @param all_years Character or numeric vector of years used for indexing.
#'
#' @return A named list of metadata items used downstream.
#' @export
#' @examples
#' \dontrun{
#' meta <- build_jags_inputs(df_idx)
#' }
build_jags_inputs <- function(df, all_years = seq(1990, 2030.5, by = 0.5)) {
  stopifnot(is.data.frame(df))
  
  # unique country-subnat mapping and counts
  index_country_subnat_tbl <- df %>%
    dplyr::group_by(Country, Region) %>% 
    dplyr::select(Country, Region, index_country, index_subnat) %>%
    dplyr::distinct() %>%
    dplyr::arrange(index_country, index_subnat)
  
  count_provinces <- index_country_subnat_tbl %>%
    dplyr::count(index_country, name = "n_subnats")
  

  # match vectors (indexing observations to parameter arrays)
  match_country <- df$index_country
  match_year <- df$index_year
  match_method <- df$index_method
  match_subnat <- df$index_subnat
  n_years <- length(all_years)
  
  # package return
  list(
    index_country_subnat_tbl = index_country_subnat_tbl,
    count_provinces = count_provinces,
    n_country = unique(index_country_subnat_tbl$Country),
    n_subnat = index_country_subnat_tbl$Region,
    match_country = match_country,
    match_year = match_year,
    match_method = match_method,
    match_subnat = match_subnat,
    n_years = n_years,
    n_obs = nrow(df)
  )
}
