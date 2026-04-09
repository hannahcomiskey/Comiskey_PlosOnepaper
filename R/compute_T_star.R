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
#' # Example: Compute T* (last observed year per subnational area) 
#' library(dplyr)
#'
#' # Create example dataset
#' df_idx <- data.frame(
#'   Country        = c("A", "A", "A", "A", "B", "B"),
#'   Region         = c("R1", "R1", "R2", "R2", "R1", "R1"),
#'   index_country  = c(1, 1, 1, 1, 2, 2),
#'   index_subnat   = c(1, 1, 2, 2, 3, 3),
#'   average_year   = c(2000, 2005, 2001, 2003, 1999, 2004),
#'   index_year     = c(1, 2, 1, 2, 1, 2)
#' )
#'
#' # Compute T*
#' T_star <- compute_T_star(df_idx)
#'
compute_T_star <- function(df) {
  required <- c("Country", "Region", "index_country", "index_subnat", "average_year", "index_year")
  stopifnot(all(required %in% names(df)))
  
  Tstar_df <- df %>%
    dplyr::group_by(Country, Region) %>%
    dplyr::filter(index_year == max(index_year, na.rm = TRUE)) %>%
    dplyr::arrange(index_subnat) %>%
    dplyr::ungroup() %>%
    dplyr::select(index_country, index_subnat, average_year, index_year) %>%
    dplyr::distinct()
  
  return(Tstar_df)
}
