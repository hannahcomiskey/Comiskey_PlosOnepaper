#' Replace Small or Missing Standard Errors
#'
#' Internal helper function that processes standard errors (SEs) in FP
#' source data. Handles zero, missing, and extremely small SE values
#' using DEFT-adjusted calculations.
#'
#' @param df A data frame containing Public.SE, Private.SE and sample sizes.
#' @param deft_lookup A DEFT lookup table.
#' @param min_n Minimum sample size threshold.
#'
#' @return A data frame with corrected SE values.
#'
#' @keywords internal
replace_small_or_missing_se <- function(df, deft_lookup, min_n = 20) {
  
  df$count_SE.NA <- is.na(df$Public.SE)
  
  normal <- df %>% dplyr::filter(count_SE.NA == 0 & Public.SE > 0.005)
  special <- df %>% dplyr::filter(Public.SE <= 0.005)
  
  if (nrow(special) > 0) {
    special <- special %>%
      dplyr::filter(Public_n >= min_n | Private_n >= min_n) %>%
      dplyr::left_join(deft_lookup, by = "average_year")
    
    for (i in seq_len(nrow(special))) {
      
      SEval <- special$Public.SE[i]
      
      if (SEval == 0) {
        # Replace 0 SE with mean of nonzero
        special$Public.SE[i] <- mean(SEval[SEval > 0], na.rm = TRUE)
      } else {
        # Compute DEFT-based replacement
        DEFT <- ifelse(is.na(special$DEFT[i]), 1.5, special$DEFT[i])
        N1 <- special$Public_n[i] + special$Private_n[i]
        phat <- 0.5 / (N1 + 1)
        SE.hat <- sqrt((phat * (1 - phat)) / N1)
        special$Public.SE[i] <- SE.hat * DEFT
      }
    }
  }
  
  df <- dplyr::bind_rows(normal, special)
  
  return(df)
}
