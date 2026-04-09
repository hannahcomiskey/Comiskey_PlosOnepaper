#' Clean Subnational Family Planning Source Data
#'
#' This function processes subnational family planning (FP) source
#' data by filtering FP2030 countries, standardising region names,
#' removing small-sample observations, fixing proportions, applying
#' lemon-squeezer transformations, and cleaning standard errors (SEs)
#' using DEFT-based imputation rules.
#'
#' @param subnat_data A data frame containing raw subnational FP source data.
#' @param deft_lookup A data frame containing DEFT values with a column
#'  named `average_year`.
#' @param fp2030_countries A character vector of FP2030 countries. A default
#'  internal vector is provided.
#' @param min_n Minimum sample size for keeping observations (default = 20).
#'
#' @return A cleaned data frame with harmonised regions, adjusted
#'  proportions, cleaned SE values, and restrictions applied.
#'
#' @details
#' Steps performed:
#' \itemize{
#'  \item Filters FP2030 countries.
#'  \item Removes regions `"NA"` and rows with insufficient sample size.
#'  \item Ensures proportions sum to 1 and replaces missing sectors with 0.
#'  \item Applies the lemon-squeezer transformation to avoid 0/1 boundaries.
#'  \item Cleans SE values
#'   \itemize{
#'    \item If one SE is 0, replaces with the mean of other sectors.
#'    \item If all SEs missing or extremely small, uses DEFT-adjusted SE.
#'   }
#'  \item Fixes region naming inconsistencies for specific countries.
#' }
#'
#' @export
#' @importFrom magrittr %>%
#' @examples
#' cleaned <- clean_fp_source_data(subnat_bivar_data)
#' 
clean_fp_source_data <- function(subnat_data,
                 deft_lookup = mcmsector::DEFT_DHS_database,
                 fp2030_countries = c(
                  "Benin", "Burkina Faso", "Cameroon",
                  "Cote d'Ivoire", "Ethiopia", "Ghana",
                  "Guinea", "Kenya", "Liberia", "Madagascar",
                  "Malawi", "Mali", "Mozambique", "Myanmar",
                  "Nepal", "Niger", "Nigeria", "Pakistan",
                  "Rwanda", "Senegal", "Tanzania", "Uganda",
                  "Zimbabwe"
                 ),
                 min_n = 20) {
  # Filter FP2030 countries and standardise region names 
  df <- subnat_data %>%
    dplyr::filter(Country %in% fp2030_countries) %>%
    dplyr::mutate(Region = stringr::str_to_title(Region)) %>%
    dplyr::filter(Region != "NA")
   
   # Keep only samples with sufficient sample size 
   df <- df %>%
     dplyr::rename(Public.SE = se.Public, Private.SE = se.Private) %>%
     dplyr::select(Country, Region, Method,  average_year, Public, Private, Public.SE, Private.SE, Public_n, Private_n) %>%
     dplyr::arrange(Country) %>%
     dplyr::filter(Public_n >= min_n | Private_n >= min_n)
   
   # Compute totals and fill missing proportions 
   df <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(check_total = sum(Public, Private, na.rm = TRUE)) %>%
    dplyr::ungroup()
   
   fill_indices <- which(df$check_total > 0.99)
   
   for (i in fill_indices) {
    na_cols <- which(is.na(df[i, c("Public", "Private")]))
    if (length(na_cols) > 0) {
     df[i, c("Public", "Private")[na_cols]] <- 0
    }
   }
   
   # Apply lemon-squeezer transformation to move away from exact 0,1 values 
   n_obs <- nrow(df)
   df <- df %>%
    dplyr::mutate(
     Public = (Public * (n_obs - 1) + 0.5) / n_obs,
     Private = (Private * (n_obs - 1) + 0.5) / n_obs
    )
   
   # Clean standard errors using helper function 
   df <- replace_small_or_missing_se(df, deft_lookup, min_n)
   
   # Fix naming inconsistencies 
   df <- fix_region_names(df)
   
   df %>%
    dplyr::arrange(Country, Region, Method, average_year)
   
   return(df)
   
}
