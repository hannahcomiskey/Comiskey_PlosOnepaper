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
  
  normal <- df %>% dplyr::filter(count_SE.NA == FALSE & Public.SE > 0.005)
  special <- df %>% dplyr::filter(count_SE.NA==TRUE | Public.SE<=0.005)
  
  if (nrow(special) > 0) {
    special <- special %>%
      dplyr::filter(Public_n >= min_n | Private_n >= min_n) %>%
      dplyr::left_join(deft_lookup %>% dplyr::rename(average_year = Year))
    
    if(nrow(special)>0) {
      col_index <- which(colnames(special)=="Public.SE")-1 # column index before SE col
      for(i in 1:nrow(special)) {
        num.SE0 <- which(special[i,c("Public.SE")]==0)
        num.SEnon0 <- which(special[i,c("Public.SE")]!=0)
        num.SEverytiny <- which(special[i,c("Public.SE")]>0 & special[i,c("Public.SE")]<0.001)
        num.SEna <- which(is.na(special[i,c("Public.SE")])==TRUE)
        DEFT <- ifelse(is.na(special$DEFT[i])==TRUE, 1.5, special$DEFT[i])
        N1 <- sum(special[i, c('Public_n', 'Private_n')], na.rm=TRUE) # Number of women surveyed
        phat <- 0.5/(N1+1) # Posterior mean of p under Jefferys prior for true prevalence of 0s.
        SE.hat <- sqrt((phat*(1-phat))/N1)
        special[i,c(col_index+num.SEna, col_index+num.SE0, col_index+num.SEverytiny)] <- SE.hat*DEFT
      }
    }
  }
  
  df <- dplyr::bind_rows(normal, special)
  
  return(df)
}
