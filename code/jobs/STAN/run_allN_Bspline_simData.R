library(rstan)
library(tidyverse)
library(tidybayes)
source('code/load_functions.R')
source('code/stan_utility.R')


# Source simulated data --------------------------------------
options(mc.cores = parallel::detectCores())

load("data/simulated_data/simulated_data_all_N_kstar_Kenya_20years.RData")

# Get logit of parameters and variance -----------------------------------------
mydata <- P_sim_df_sample[,c("Public")] %>%
  rowwise() %>%
  mutate(Public = ifelse(Public < 0.0001 , 0, Public))

logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)))

# # testing splines ------------------------------------------------------------
all_years <- seq(1:max(P_sim_df_sample$index_year))
B <- splines::bs(all_years, degree=3, knots=c(seq(1, 50, by=4)),  intercept = TRUE) # creating the B-splines
num_data <- length(all_years)
num_basis <- nrow(B)

# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- rep(1, 6)

n_all_years <- length(all_years)
M_count = 5

# Run stan model ------------------

sm<-stan_model("model/STAN/fit_basis.stan")

fit<-sampling(sm,
              iter=5000,
              warmup=1000,
              thin=2,
              chains=3,
              data = list(Y = as.vector(unlist(logit.data[,c("logit.Public")])),
                          X = all_years,
                          num_years = n_all_years,
                          num_data = nrow(logit.data),
                          num_knots = ncol(B),
                          P_count = P,
                          M_count = M,
                          S_count = 2,
                          B = B,
                          matchsubnat = simmatchsubnat,
                          matchmethod = simmatchmethod,
                          matchyears = simmatchyears), 
              control = list(adapt_delta=0.95))

saveRDS(fit, file='results/STAN/STAN_model_fit_basis_SimKenya20.RDS')
