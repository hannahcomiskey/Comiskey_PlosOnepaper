#' Standardise country names and attach super-region classification
#'
#' @description
#' Harmonises a few common country name variants and joins a country-area
#' classification table to produce a `Super_region` column.
#'
#' @param df A data frame with a `Country` column.
#' @param area_classification A data frame with columns `Country or area` and `Region`
#'   (e.g., `mcmsupply::Country_and_area_classification`).
#' @param country_name_fixes Named character vector of replacements (optional).
#'
#' @return The input data frame with `Super_region` added and country names fixed.
#' @export
#' @examples
#' \dontrun{
#' ac <- mcmsupply::Country_and_area_classification
#' df2 <- standardize_country_names(FP_source_data_wide, ac)
#' }
standardize_country_names <- function(df,
                                      area_classification,
                                      country_name_fixes = c(
                                        "Bolivia (Plurinational State of)" = "Bolivia",
                                        "Republic of Moldova" = "Moldova",
                                        "Viet Nam" = "Vietnam",
                                        "Kyrgyzstan" = "Kyrgyz Republic"
                                      )) {
  stopifnot(is.data.frame(df), is.data.frame(area_classification))
  df <- df %>%
    dplyr::mutate(Country = as.character(Country))
  
  # apply name fixes if present
  if (length(country_name_fixes) > 0) {
    df$Country <- dplyr::recode(df$Country, !!!country_name_fixes)
  }
  
  # prepare classification table
  ac <- area_classification %>%
    dplyr::select(`Country or area`, Region) %>%
    dplyr::rename(Country = `Country or area`, Super_region = Region) %>%
    dplyr::mutate(Country = as.character(Country))
  
  # join and return
  df %>%
    dplyr::left_join(ac, by = "Country")
  
  return(df)
}
