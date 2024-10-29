library(rstan)
library(tidyverse)
library(tidybayes)


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

# Source simulated data --------------------------------------
options(mc.cores = parallel::detectCores())

load("data/simulated_data/simulated_data_all_N_kstar_Kenya.RData")

# Get logit of parameters and variance -----------------------------------------
mydata <- P_sim_df_sample[,c("Public")] %>%
  mutate(Public.SE = 0.1) %>%
  rowwise() %>%
  mutate(Public = ifelse(Public < 0.0001 , 0, Public))

logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))

# # testing splines ------------------------------------------------------------
all_years <- -10:30
B <- bs_bbase_precise(all_years)
Bik <- B$B.ik
K <-dim(Bik)[2]
H = K-1
kstar = B$Kstar

year_index_table <- tibble(Year = all_years, index_year = 1:length(all_years))

P_sim_df_sample <- P_sim_df_sample %>% 
  rename(Year = index_year) %>%
  mutate_if(is.character, as.numeric) %>%
  left_join(year_index_table)

# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- matchcountry
n_all_years <- length(all_years)
M_count = 5

Bik_array <- array(Bik, dim=c(P,n_all_years,K))

## The required data ------------------------------
inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  Bik = Bik_array,
                  n_years = n_all_years,
                  n_obs = nrow(logit.data),
                  K = K,
                  H = H,
                  kstar = rep(kstar, P),
                  zero = 0,
                  C_count = C,
                  P_count = P,
                  M_count = M,
                  S_count = 2,
                  matchsubnat = simmatchsubnat,
                  matchcountry = simmatchcountry,
                  matchmethod = simmatchmethod,
                  matchyears = simmatchyears)

## Parameters to look at ------------------------------
pars <- c("alpha_pms", # required for P
          "delta_k",
          "beta_c",
          "beta_k",
          "sigma_alpha",
          "sigma_beta",
          "sigma_delta",
          "P")

# Run stan model ------------------

fit <- stan(
  file = "model/STAN/kstar_spline_public_private_logit_allN.stan",  # Stan program
  data = inputdata,    # named list of data
  pars = pars,
  iter = 40000,         # total number of iterations per chain
  warmup = 10000,
  thin=15,
  chains=3,
  save_warmup = FALSE,
  control=list(adapt_delta=0.99)
)


saveRDS(fit, file='results/STAN_model_test_kstar_simdata_N_alpha_N_delta_KenyaSim.RDS')
