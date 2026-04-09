#' Harmonise Region Names for Selected Countries
#'
#' This function standardises region names for Burkina Faso, Rwanda,
#' Nigeria, and Cote d'Ivoire to ensure consistent naming across years.
#'
#' @param df A data frame containing `Country` and `Region` columns.
#'
#' @return A data frame with modified region names.
#'
#' @export
#'
#' @examples
#' clean <- fix_region_names(mcmsector::subnat_bivar_data)
#' 
fix_region_names <- function(df) {
  
  fix_country <- function(data, country, replacements) {
    tmp <- data %>%
      dplyr::filter(Country == country) %>%
      dplyr::mutate(Region = dplyr::case_when(!!!replacements))
    rest <- data %>% dplyr::filter(Country != country)
    dplyr::bind_rows(rest, tmp)
  }
  
  # Burkina Faso
  df <- fix_country(df, "Burkina Faso", rlang::exprs(
    Region == "North" ~ "Nord",
    Region == "East" ~ "Est",
    Region == "West" ~ "Centre-Ouest",
    Region == "Central/South" ~ "Centre-Sud",
    Region == "Ouagadougou" ~ "Centre Including Ouagadougou",
    TRUE ~ Region
  ))
  
  # Rwanda
  df <- fix_country(df, "Rwanda", rlang::exprs(
    Region %in% c("Ouest") ~ "West",
    Region == "Nord" ~ "North",
    Region == "Sud" ~ "South",
    Region == "Est" ~ "East",
    Region %in% c("Ville de Kigali", "Ville De Kigali", "Kigali City") ~ "Kigali",
    Region == "Butare, Gitarama (Central, South)" ~ "South",
    Region == "Cyangugu, Gikongoro (Southwest)" ~ "South",
    Region == "Byumba, Kibungo, Umutara (Northeast)" ~ "North",
    Region == "Gisenyi, Kibuye, Ruhengeri (Northwest)" ~ "West",
    TRUE ~ Region
  ))
  
  # Nigeria
  df <- fix_country(df, "Nigeria", rlang::exprs(
    Region == "Northeast" ~ "North East",
    Region == "Northwest" ~ "North West",
    Region == "Southeast" ~ "South East",
    Region == "Southwest" ~ "South West",
    Region == "Central" ~ "North Central",
    TRUE ~ Region
  ))
  
  # Cote d'Ivoire
  df <- fix_country(df, "Cote d'Ivoire", rlang::exprs(
    Region == "Center-East" ~ "Center East",
    Region == "Center-North" ~ "Center North",
    Region == "Center-West" ~ "Center West",
    Region == "Center-South" ~ "Center South",
    TRUE ~ Region
  ))
  
  return(df)
}
