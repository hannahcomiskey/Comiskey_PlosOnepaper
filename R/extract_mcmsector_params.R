#' Extract alpha_pms and beta.k MCMC samples from a JAGS model object
#'
#' @description
#' Given a JAGS model object containing `BUGSoutput$sims.array`, this function
#' extracts the posterior simulations for `alpha_pms` and `beta.k` into
#' clean matrices suitable for downstream computation.
#'
#' @param mod A JAGS model object loaded with `readRDS()`, containing
#'   `mod$BUGSoutput$sims.array`.
#' @param n_method Integer; number of contraceptive methods.
#' @param n_subnat Integer; number of subnational regions.
#' @param n_beta Integer; number of beta coefficients (typically 13).
#'
#' @return A list with:
#' \describe{
#'   \item{alpha_pms}{A matrix of posterior samples for alpha parameters.}
#'   \item{beta_k}{A matrix of posterior samples for beta parameters.}
#' }
#' @export
#'
#' @examples
#' #  Example: Extract parameters from mock JAGS output 
#' set.seed(123)
#'
#' # Create a fake sims.array
#' n_iter   <- 10
#' n_chain  <- 2
#' n_method <- 2
#' n_subnat <- 3
#' n_beta   <- 2
#'
#' var_names <- c(
#'   paste0("alpha_pms[", rep(1:n_method, each = n_subnat),
#'          ",", rep(1:n_subnat, times = n_method), "]"),
#'   paste0("beta.k[",
#'          rep(1:n_method, each = n_subnat * n_beta), ",",
#'          rep(rep(1:n_subnat, each = n_beta), times = n_method), ",",
#'          rep(1:n_beta, times = n_method * n_subnat), "]")
#' )
#'
#' sims_array <- array(
#'   rnorm(n_iter * n_chain * length(var_names)),
#'   dim = c(n_iter, n_chain, length(var_names)),
#'   dimnames = list(NULL, NULL, var_names)
#' )
#'
#' # Mock JAGS object
#' mod <- list(BUGSoutput = list(sims.array = sims_array))
#'
#' # Run function
#' params <- extract_mcmsector_params(
#'   mod,
#'   n_method = n_method,
#'   n_subnat = n_subnat,
#'   n_beta   = n_beta
#' )
#' 
extract_mcmsector_params <- function(mod, n_method, n_subnat, n_beta = 13) {
  
  sims <- mod$BUGSoutput$sims.array
  vars <- dimnames(sims)[[3]]
  
  # alpha --
  alpha_start <- grep("alpha_pms\\[1,1\\]", vars)[1]
  alpha_end   <- grep(sprintf("alpha_pms\\[%s,%s\\]", n_method, n_subnat), vars)[1]
  
  alpha_raw <- sims[ , , alpha_start:alpha_end ]
  alpha_mat <- rbind(alpha_raw[, 1, ], alpha_raw[, 2, ])
  
  # beta ---
  beta_start <- grep("beta.k\\[1,1,1\\]", vars)[1]
  beta_end   <- grep(sprintf("beta.k\\[%s,%s,%s\\]", n_method, n_subnat, n_beta), vars)[1]
  
  beta_raw <- sims[ , , beta_start:beta_end ]
  beta_mat <- rbind(beta_raw[, 1, ], beta_raw[, 2, ])
  
  return(list(
    alpha_pms = alpha_mat,
    beta_k = beta_mat
  ))
}
