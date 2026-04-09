#' Build a B-spline basis for each subnational area
#'
#' @description
#' For each subnational area this builds a b-spline basis matrix across `all_years`
#' using `bs_bbase_precise()` which is a small wrapper around `splines::bs()`.
#' The returned arrays mimic the shape used in your original script.
#'
#' @param T_star Data frame returned by `compute_T_star()` with `average_year` and `index_subnat`.
#' @param all_years Numeric vector of years for which the basis is evaluated.
#' @param nseg Integer number of segments for interior knots (default 10).
#'
#' @return A list with elements:
#' \itemize{
#'  \item B_ik: array (n_subnat x length(all_years) x K)
#'  \item Kstar: integer vector length n_subnat (effective K for each area)
#'  \item knots_all: matrix of knots (n_subnat x K)
#'  \item K: number of basis columns
#'  \item H: K - 1
#' }
#' @export
#' @examples
#' # Example: Build B-spline basis 
#' library(dplyr)
#'
#' # Create a small T_star dataset (as returned by compute_T_star)
#' T_star <- data.frame(
#'   index_country = c(1, 1, 2),
#'   index_subnat  = c(1, 2, 3),
#'   average_year  = c(2005, 2010, 2008),
#'   index_year    = c(2, 3, 2)
#' )
#'
#' # Define years to evaluate the spline basis
#' all_years <- seq(2000, 2015, by = 1)
#'
#' # Build spline basis
#' splines_out <- build_spline_basis(
#'   T_star    = T_star,
#'   all_years = all_years,
#'   nseg      = 5
#' )
#' 
build_spline_basis <- function(T_star, all_years = seq(1990, 2030.5, by = 0.5), nseg = 10) {
  stopifnot(is.data.frame(T_star))
  n_subnat <- nrow(T_star)
  K <- nseg + 3  # cubic B-spline default (nseg interior knots + degree)
  B_ik <- array(0, dim = c(n_subnat, length(all_years), K))
  knots_all <- matrix(NA_real_, nrow = n_subnat, ncol = nseg + 3)
  Kstar <- integer(n_subnat)
  
  for (i in seq_len(n_subnat)) {
    lastobs <- T_star$average_year[i]
    res <- bs_bbase_precise(x = all_years, lastobs = lastobs, nseg = nseg)
    # res$B is matrix length(all_years) x K
    B_ik[i, , ] <- as.matrix(res$B)
    Kstar[i] <- res$Kstar
    knots_all[i, seq_len(length(res$knots.k))] <- res$knots.k
  }
  
  list(
    B_ik = B_ik,
    Kstar = Kstar,
    knots_all = knots_all,
    K = K,
    H = K - 1
  )
}

#' Small wrapper to create a B-spline basis similar to the original bs_bbase_precise
#'
#' @param x numeric vector of times to evaluate the basis
#' @param lastobs numeric last observed year (upper knot placement)
#' @param nseg integer number of interior segments (default 10)
#'
#' @return list(B.ik = basis matrix (length(all_years) x K), Kstar = effective K, knots.k = knots)
#' @keywords internal
bs_bbase_precise <- function(x = x,lastobs = max(x), xl = min(x), xr = max(x), nseg = 10, deg = 3) {
  # Compute the length of the partitions
  dx <- (xr - xl) / nseg
  # Compute position of knot before last observation
  dk <- lastobs
  # Create equally spaced knots
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  # Find index of closest knot to dk
  dk_index <- which.min(abs(knots-dk))
  # Find transformation to knot placement so that dk is a knot 
  ktrans <- (dk-knots)[dk_index]
  # Add transformation to knots
  knotsnew <- knots + ktrans
  # Use bs() function to generate the B-spline basis
  get_bs_matrix <- matrix(splines::bs(x, knots = knotsnew, degree = deg, Boundary.knots = c(knotsnew[1], knotsnew[length(knotsnew)])), nrow = length(x))
  
  # Remove columns that contain zero only
  bs_matrix <- get_bs_matrix[, -c(1:deg, ncol(get_bs_matrix):(ncol(get_bs_matrix) - deg))]
  
  used_knots <- knotsnew[-c(1,2,length(knotsnew),(length(knotsnew)-1))]
  Kstar <- which(used_knots==dk)
  
  return(list(B.ik = bs_matrix, ##<< Matrix, each row is one observation, each column is one B-spline.
              knots.k = used_knots, ##<< Vector of transformed knots.
              Kstar = Kstar # Knot point of last observation
  ))
}

