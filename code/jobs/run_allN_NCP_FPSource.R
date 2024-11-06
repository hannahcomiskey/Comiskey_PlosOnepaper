library(rstan)
library(tidyverse)
library(tidybayes)

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

# Source simulated data --------------------------------------
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## The required data ------------------------------
inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                  Bik = B.ik,
                  n_years = n_all_years,
                  n_obs = nrow(logit.data),
                  K = K,
                  H = H,
                  kstar = Kstar,
                  zero = 0,
                  C_count = length(n_country),
                  P_count = length(n_subnat),
                  M_count = length(n_method),
                  S_count = 2,
                  matchsubnat = FP_source_data_wide$index_subnat,
                  matchcountry = index_country_subnat_tbl$index_country,
                  matchmethod = FP_source_data_wide$index_method,
                  matchyears = FP_source_data_wide$index_year)

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
  file = "model/STAN/ncp_allN_model_SE.stan",  # Stan program
  data = inputdata,    # named list of data
  pars = pars,
  iter = 10000,         # total number of iterations per chain
  warmup = 2000,
  thin=4,
  chains=3,
  save_warmup = FALSE,
  control=list(adapt_delta=0.99)
)


saveRDS(fit, file='results/STAN_model_Kenya_allN_NCP_kstar.stan')
