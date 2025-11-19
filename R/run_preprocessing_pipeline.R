#' High-level preprocessing pipeline
#'
#' @description
#' Master wrapper that runs the full preprocessing pipeline:
#'  - standardise country names & join area classification
#'  - add index variables
#'  - build JAGS metadata
#'  - compute T_star
#'  - build spline bases
#'
#' This returns a list containing cleaned data, metadata and spline objects.
#'
#' @param raw_df Raw FP source data frame (unindexed).
#' @param area_classification A country/area classification table (e.g. mcmsupply data).
#' @param deft_lookup Optional DEFT lookup (passed through to SE cleaning if used earlier).
#' @param methods Character vector of methods.
#' @param all_years Numeric vector of years to index across.
#' @param nseg Number of spline segments.
#'
#' @return A named list with elements:
#' \itemize{
#'  \item data: cleaned + indexed data frame
#'  \item meta: jags metadata list
#'  \item T_star: final-observation table
#'  \item splines: spline basis object
#' }
#' @export
#' @examples
#' \dontrun{
#' out <- run_preprocessing_pipeline(raw_df = FP_source_data_wide,
#'                                   area_classification = mcmsupply::Country_and_area_classification)
#' }
run_preprocessing_pipeline <- function(raw_df,
                                       area_classification,
                                       deft_lookup = NULL,
                                       methods = c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills"),
                                       all_years = seq(1990, 2030.5, by = 0.5),
                                       nseg = 10,
                                       min_n = 20) {
  # 1. standardise country names & attach super-region
  df1 <- standardize_country_names(raw_df, area_classification)
  
  # 2. add index variables
  df2 <- add_index_variables(df1, methods = methods, all_years = all_years)
  
  # Optional: if user wants SE cleaning before indexing, they can run the earlier cleaning functions
  # (we kept this pipeline focused on indexing & metadata)
  # 3. Build jags metadata
  meta <- build_jags_metadata(df2, all_years = all_years)
  
  # 4. Compute T_star
  T_star <- compute_T_star(df2)
  
  # 5. Build splines
  splines <- build_spline_basis(T_star, all_years = all_years, nseg = nseg)
  
  return(list(
    data = df2,
    meta = meta,
    T_star = T_star,
    splines = splines
  ))
}
