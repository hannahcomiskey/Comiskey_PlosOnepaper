#' Create index variables for subnational modelling
#'
#' @description
#' Adds integer index columns for subnational area (`index_subnat`),
#' country (`index_country`), method (`index_method`) and year
#' (`index_year`) to facilitate modeling and JAGS data construction.
#'
#' @param df A data frame that contains at least `Country`, `Region`,
#'   `Method` and `average_year` columns.
#' @param methods Character vector of methods in the desired order (default examples given).
#' @param all_years Numeric vector of years to index against (defaults to 1990–2030 by 0.5).
#'
#' @return The input data frame augmented with index columns.
#' @export
#' @examples
#' \dontrun{
#' df_idx <- add_index_variables(FP_source_data_wide)
#' }
add_index_variables <- function(df,
                                methods = c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills"),
                                all_years = seq(1990, 2030.5, by = 0.5)) {
  stopifnot(is.data.frame(df))
  required <- c("Country", "Region", "Method", "average_year")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols)) {
    stop("Input data frame must contain columns: ", paste(required, collapse = ", "))
  }
  
  # create country - subnat unique table
  country_subnat_tbl <- df %>%
    dplyr::select(Country, Region) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Country, Region) %>%
    dplyr::mutate(index_country = as.integer(as.factor(Country)),
                  index_subnat = dplyr::row_number())
  
  # add indices to main df using left_join to keep row order
  df2 <- df %>%
    dplyr::left_join(country_subnat_tbl, by = c("Country", "Region"))
  
  # country index (factor level integer)
  df2 <- df2 %>%
    dplyr::mutate(index_country = as.integer(factor(Country, levels = unique(country_subnat_tbl$Country))),
                  index_method = as.integer(factor(Method, levels = methods)),
                  index_year = match(average_year, all_years))
  
  # If any years are not matched, warn and create NA => user must ensure all_years covers average_year
  if (any(is.na(df2$index_year))) {
    warning("Some average_year values were not found in 'all_years'. Please expand all_years or check values.")
  }
  
  return(df2)
}
