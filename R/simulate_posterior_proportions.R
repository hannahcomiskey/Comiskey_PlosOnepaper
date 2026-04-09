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
#'   \item{P}{A 5D array of public/private probabilities.}
#' }
#' @examples
#' set.seed(123)
#'
#' # Dimensions
#' n_samps  <- 100
#' n_method <- 2
#' n_subnat <- 3
#' n_years  <- 5
#' n_beta   <- 4
#'
#' all_years <- 2000:(2000 + n_years - 1)
#'
#' # Simulate alpha posterior samples
#' alpha_pms <- matrix(rnorm(n_samps * n_method * n_subnat),
#'                     nrow = n_samps)
#'
#' colnames(alpha_pms) <- as.vector(
#'   outer(
#'     paste0("alpha_pms[", 1:n_method, ","),
#'     paste0(1:n_subnat, "]"),
#'     paste0
#'   )
#' )
#'
#' # Simulate beta posterior samples 
#' beta_k <- matrix(rnorm(n_samps * n_method * n_subnat * n_beta),
#'                  nrow = n_samps)
#'
#' colnames(beta_k) <- unlist(
#'   lapply(1:n_method, function(m) {
#'     lapply(1:n_subnat, function(p) {
#'       paste0("beta.k[", m, ",", p, ",", 1:n_beta, "]")
#'     })
#'   })
#' )
#'
#' # Create basis function array 
#' B_ik <- array(runif(n_subnat * n_years * n_beta),
#'               dim = c(n_subnat, n_years, n_beta))
#'
#' # Run simulation 
#' P <- simulate_posterior_proportions(
#'   alpha_pms = alpha_pms,
#'   beta_k    = beta_k,
#'   B_ik      = B_ik,
#'   n_method  = n_method,
#'   n_subnat  = n_subnat,
#'   all_years = all_years,
#'   n_samps   = n_samps
#' )
#'
#' @export
simulate_posterior_proportions <- function(alpha_pms,
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
      
      alpha_samp <- alpha_pms[1:n_samps, grepl(sprintf("alpha_pms\\[%s,%s\\]", m, p),
                                       colnames(alpha_pms)) ]
      
      beta_samp <- beta_k[ , grepl(sprintf("beta.k\\[%s,%s,", m, p),
                                   colnames(beta_k)) ]
      
      for (t in seq_len(n_years)) {
        xb <- rowSums(beta_samp[1:n_samps, 1:n_beta, drop=FALSE] *
                        matrix(B_ik[p, t, ], nrow=n_samps, ncol=n_beta, byrow=TRUE))
        
        z[ , m, p, t] <- alpha_samp + xb
        
        P[ , 1, m, p, t] <- 1 / (1 + exp(-z[ , m, p, t]))
        P[ , 2, m, p, t] <- 1 - P[ , 1, m, p, t]
      }
    }
  }
  
  return(P)
}
