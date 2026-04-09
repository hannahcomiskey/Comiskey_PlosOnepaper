#' Summarise posterior probability samples into mean and quantiles
#' @importFrom dplyr across everything
#' @description
#' Converts a 5D posterior draws array `P` into a tidy long data frame containing
#' posterior means and 95% credible intervals for each sector, method, region,
#' and year.
#' @importFrom stats quantile
#' @param P A 5D array of posterior samples with dimensions:
#'   (sample, sector, method, subnat, year).
#' @param method_index_table Tibble mapping method indices to method names.
#' @param sector_index_table Tibble mapping sector indices to sector names.
#' @param subnat_index_table Mapping subnational indices to regions/countries.
#' @param year_index_table Mapping year index to calendar year.
#'
#' @return A tibble with posterior mean, lower 95%, and upper 95% intervals.
#' @export
#' @examples
#' # Example: Summarise posterior probabilities
#' library(dplyr)
#' library(tidyr)
#'
#' # Simulate a small 5D posterior array:
#' # Dimensions: sample x sector x method x subnat x year
#' n_samps  <- 10
#' n_sector <- 2
#' n_method <- 2
#' n_subnat <- 3
#' n_year   <- 4
#'
#' set.seed(123)
#' P <- array(runif(n_samps * n_sector * n_method * n_subnat * n_year),
#'            dim = c(n_samps, n_sector, n_method, n_subnat, n_year))
#'
#' # Create minimal index tables
#' method_tbl <- tibble(index_method = 1:n_method, Method = c("MethodA", "MethodB"))
#' sector_tbl <- tibble(index_sector = 1:n_sector, SectorName = c("Public", "Private"))
#' subnat_tbl <- tibble(index_subnat = 1:n_subnat, Region = paste0("Region", 1:n_subnat))
#' year_tbl   <- tibble(index_year = 1:n_year, Year = 2000 + 0:(n_year-1))
#'
#' # Summarise posterior
#' df_summary <- summarise_posterior_proportions(
#'   P,
#'   method_index_table = method_tbl,
#'   sector_index_table = sector_tbl,
#'   subnat_index_table = subnat_tbl,
#'   year_index_table   = year_tbl
#' )
#' 
summarise_posterior_proportions <- function( P,
                                             method_index_table,
                                             sector_index_table,
                                             subnat_index_table,
                                             year_index_table) {
  
  # Means -----------------------------------------------------------
  P_mean <- apply(P, c(2,3,4,5), mean)
  P_mean <- plyr::adply(P_mean, .margins = c(2, 3, 4))
  colnames(P_mean) <- c("index_method", "index_subnat", "index_year", "Public", "Private")
  
  P_mean <- P_mean %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    tidyr::pivot_longer(cols = c("Public", "Private"), names_to="Sector", values_to="Mean") %>%
    dplyr::left_join(method_index_table, by="index_method")
  
  # Quantiles --------------------------------------------------------
  P_q <- apply(P, c(2,3,4,5), quantile, probs=c(0.025, 0.975))
  P_q <- plyr::adply(P_q, .margins = c(2,3,4,5))
  colnames(P_q) <- c("index_sector", "index_method", "index_subnat", "index_year",
                     "lower_95", "upper_95")
  
  P_q <- P_q %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    dplyr::left_join(method_index_table, by="index_method") %>%
    dplyr::left_join(sector_index_table, by="index_sector") %>%
    dplyr::left_join(subnat_index_table, by=c("index_subnat")) %>%
    dplyr::left_join(year_index_table, by="index_year")
  
  # Combine ---------------------------------------------------------
  P_df <- dplyr::left_join(P_mean, P_q)
  return(P_df)
}
