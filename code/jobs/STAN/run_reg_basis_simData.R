library(rstan)
library(tidyverse)
library(tidybayes)
library(splines)
source('code/load_functions.R')

# Source simulated data --------------------------------------------------------
options(mc.cores = parallel::detectCores())

load("data/simulated_data/simulated_data_all_N_kstar_Kenya_50years.RData")

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

Y = logit.data$logit.Public

## testing splines -------------------------------------------------------------
X <- seq(from=-5, to=55, by=1) # generating inputs

B <- t(bs(X, knots=seq(-5,55,5), degree=3, intercept = TRUE)) # creating the B-splines

# B <- bs_bbase_precise(X)
# B <- t(B$B.ik)

num_years <- length(X)
num_basis <- nrow(B)
num_data <- length(Y)

year_index_table <- tibble(Year = X, index_year = 1:length(X))

P_sim_df_sample <- P_sim_df_sample %>% 
  rename(Year = index_year) %>%
  mutate_if(is.character, as.numeric) %>%
  left_join(year_index_table)

## Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- matchcountry
n_all_years <- length(all_years)
M_count = 5

## The required data -----------------------------------------------------------
sm<-stan_model("model/STAN/fit_reg_basis.stan")

fit<-sampling(sm,
              iter=5000, 
              data = list(num_data = num_data, 
                          num_basis = num_basis, 
                          num_years = num_years,
                          Y = Y, 
                          X=X, 
                          B=B,
                          P_count = P,
                          M_count = M,
                          S_count = 2,
                          matchsubnat = simmatchsubnat,
                          matchcountry = simmatchcountry,
                          matchmethod = simmatchmethod,
                          matchyears = simmatchyears,
                          v = 6,
                          s_sq = 10), 
              control = list(adapt_delta=0.95))

saveRDS(fit, file='results/STAN_model_reg_basis_KenyaSim_50years_test.RDS')


