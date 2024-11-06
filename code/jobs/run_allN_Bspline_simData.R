library(rstan)
library(tidyverse)
library(tidybayes)
source('code/load_functions.R')

# Source simulated data --------------------------------------
options(mc.cores = parallel::detectCores())

load("data/simulated_data/simulated_data_all_N_kstar_Kenya_new.RData")

# Get logit of parameters and variance -----------------------------------------
mydata <- P_sim_df_sample[,c("Public")] %>%
  rowwise() %>%
  mutate(Public = ifelse(Public < 0.0001 , 0, Public))

logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)))

# # testing splines ------------------------------------------------------------
all_years <- seq(1:max(P_sim_df_sample$index_year))
B <- t(splines::bs(all_years, degree=3, knots=c(seq(1, 20, by=4)),  intercept = TRUE)) # creating the B-splines
num_data <- length(all_years)
num_basis <- nrow(B)

# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- rep(1, 6)

n_all_years <- length(all_years)
M_count = 5

# let p0=3, D=H and n=n_obs 
# See 3.12 Vehtari paper.
scale_global = 3/((H-3)*sqrt(nrow(logit.data)))

## The required data ------------------------------
inputdata <- list(Y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  X = all_years,
                  num_years = n_all_years,
                  num_data = nrow(logit.data),
                  num_knots = nrow(B),
                  P_count = P,
                  M_count = M,
                  S_count = 2,
                  B = B,
                  matchsubnat = simmatchsubnat,
                  matchmethod = simmatchmethod,
                  matchyears = simmatchyears,
                  scale_global=scale_global
                  )

## Parameters to look at ------------------------------
pars <- c("a",
          "lambda",
          "beta_c",
          "a0",
          "sigma_a0",
          "P")

# Run stan model ------------------

fit <- stan(
  file = "model/STAN/regularised_Bspline_model.stan",  # Stan program
  data = inputdata,    # named list of data
  pars = pars,
  iter = 10000,         # total number of iterations per chain
  warmup = 2000,
  thin=4,
  chains=3,
  save_warmup = FALSE,
  control=list(adapt_delta=0.99)
)


saveRDS(fit, file='results/STAN_model_Bspline_allN_SimKenya.RDS')
