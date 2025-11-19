#' Compute T* (last observed year per subnational area)
#'
#' @description
#' For each country-region (subnational) pair, find the last observed index_year
#' and the associated average_year.
#'
#' @param df Data frame containing at least `Country`, `Region`, `index_country`, `index_subnat`,
#'   `average_year` and `index_year`.
#'
#' @return A data frame with one row per subnational area containing its final observation information.
#' @export
#' @examples
#' \dontrun{
#' T_star <- compute_T_star(df_idx)
#' }
compute_T_star <- function(df) {
  required <- c("Country", "Region", "index_country", "index_subnat", "average_year", "index_year")
  stopifnot(all(required %in% names(df)))
  
  df %>%
    dplyr::group_by(Country, Region) %>%
    dplyr::filter(index_year == max(index_year, na.rm = TRUE)) %>%
    dplyr::arrange(index_subnat) %>%
    dplyr::ungroup() %>%
    dplyr::select(index_country, index_subnat, average_year, index_year) %>%
    dplyr::distinct()
  
  return(df)
}
