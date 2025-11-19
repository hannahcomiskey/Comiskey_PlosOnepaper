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
  
  # year sequence for integer years (floor of index_year)
  t_seq_2 <- floor(match_year)
  year_seq <- seq(min(t_seq_2, na.rm = TRUE), max(t_seq_2, na.rm = TRUE), by = 1)
  n_years <- length(year_seq)
  
  # package return
  list(
    index_country_subnat_tbl = index_country_subnat_tbl,
    count_provinces = count_provinces,
    match_country = match_country,
    match_year = match_year,
    match_method = match_method,
    match_subnat = match_subnat,
    t_seq_2 = t_seq_2,
    year_seq = year_seq,
    n_years = n_years,
    n_obs = nrow(df)
  )
}
