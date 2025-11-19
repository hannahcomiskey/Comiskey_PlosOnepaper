#' Simulate posterior public/private proportions from alpha and beta parameters
#'
#' @description
#' Computes posterior draws of the logit-scale proportions and the corresponding
#' public/private probabilities for all methods, regions, and years.
#'
#' @param alpha_pms Matrix of alpha posterior samples.
#' @param beta_k Matrix of beta posterior samples.
#' @param B_ik A 3D array of basis functions with dimensions:
#'   (subnational, year, n_beta).
#' @param n_method Number of methods.
#' @param n_subnat Number of subnational units.
#' @param all_years Vector of years included in the model.
#' @param n_samps Number of posterior samples to compute (default = 4000).
#'
#' @return A list containing:
#' \describe{
#'   \item{z}{A 4D array of logit-scale draws.}
#'   \item{P}{A 5D array of public/private probabilities.}
#' }
#'
#' @export
simulate_posterior_proportions <- function( alpha_pms,
                                            beta_k,
                                            B_ik,
                                            n_method,
                                            n_subnat,
                                            all_years,
                                            n_samps = 4000) {
  
  n_years <- length(all_years)
  n_beta <- dim(B_ik)[3]
  
  z <- array(NA, dim = c(n_samps, n_method, n_subnat, n_years))
  P <- array(NA, dim = c(n_samps, 2, n_method, n_subnat, n_years))
  
  for (p in seq_len(n_subnat)) {
    for (m in seq_len(n_method)) {
      
      alpha_samp <- alpha_pms[ , grepl(sprintf("alpha_pms\\[%s,%s\\]", m, p),
                                       colnames(alpha_pms)) ]
      
      beta_samp <- beta_k[ , grepl(sprintf("beta.k\\[%s,%s,", m, p),
                                   colnames(beta_k)) ]
      
      for (t in seq_len(n_years)) {
        
        # B.ik[p, t, ] %*% beta_samp[s, 1:n_beta] computed fast via matrix mult
        xb <- rowSums(beta_samp[ , 1:n_beta, drop=FALSE] *
                        matrix(B_ik[p, t, ], nrow=n_samps, ncol=n_beta, byrow=TRUE))
        
        z[ , m, p, t] <- alpha_samp + xb
        
        P[ , 1, m, p, t] <- 1 / (1 + exp(-z[ , m, p, t]))
        P[ , 2, m, p, t] <- 1 - P[ , 1, m, p, t]
      }
    }
  }
  
  return(list(z = z, P = P))
}
