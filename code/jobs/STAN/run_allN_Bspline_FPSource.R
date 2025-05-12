library(rstan)
library(tidyverse)
library(tidybayes)
library(splines)
source('code/load_functions.R')

# Source code --------------------------------------
source("code/load_functions.R")
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Get logit of parameters and variance --------------------
mydata <- FP_source_data_wide[,c("Public", "Public.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))

Y = logit.data$logit.Public

## testing splines -------------------------------------------------------------

B <- t(bs(all_years, knots=seq(min(all_years),max(all_years),5), degree=3, intercept = TRUE)) # creating the B-splines

num_years <- length(all_years)
num_basis <- nrow(B)
num_data <- length(Y)

## The required data -----------------------------------------------------------
sm<-stan_model("model/STAN/fit_basis_expanded.stan")

fit<-sampling(sm,
              iter=12000, 
              warmup=4000,
              thin = 4,
              chains=3,
              data = list(num_data = num_data, 
                          num_basis = num_basis, 
                          num_years = num_years,
                          Y = Y, 
                          X=all_years, 
                          B=B,
                          P_count = length(n_subnat),
                          M_count = length(n_method),
                          S_count = 2,
                          matchsubnat = FP_source_data_wide$index_subnat,
                          matchmethod = FP_source_data_wide$index_method,
                          matchyears = FP_source_data_wide$index_year), 
              control = list(adapt_delta=0.99),
              save_warmup = FALSE
              )

saveRDS(fit, file='results/STAN_model_basis_expand_4country.RDS')
