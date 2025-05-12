library(rstan)
library(tidyverse)
library(tidybayes)

# Source simulated data --------------------------------------------------------
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load("data/simulated_data_factorcov_new.RData")

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
B <- splines::bs(all_years, df=10, degree=3, intercept = FALSE)
K <-dim(B)[2]
B.ik <- B

# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- matchcountry
n_all_years <- length(all_years)

D=2
M_count = 5
OD_count  = D*(M_count-D)+ D*(D-1)/2

## The required data -----------------------------------------------------------
inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                  Bik = B.ik,
                  intercept = rep(1, n_all_years),
                  n_years = n_all_years,
                  a_mu = rep(0, 5),
                  beta_mu = rep(0, 5),
                  n_obs = nrow(logit.data),
                  K = K,
                  D=2,
                  OD_count = OD_count,
                  C_count = C,
                  P_count = P,
                  M_count = 5,
                  S_count = 2,
                  matchsubnat = simmatchsubnat,
                  matchcountry = simmatchcountry,
                  matchmethod = simmatchmethod,
                  matchyears = simmatchyears)


## Parameters to look at -------------------------------------------------------
pars <- c("L",
          "Q",
          "L_t",
          "L_d",
          "psi",
          "a_k",
          "tau_a",
          "a_raw",
          "P",
          "beta_c",
          "alpha_pms")

# Run stan model ---------------------------------------------------------------

fit <- stan(
  file = "model/STAN/Bspline_factorcov_public_private_logit.stan",  # Stan program
  data = inputdata,    # named list of data
  pars = pars,
  iter = 40000,         # total number of iterations per chain
  warmup = 10000,
  thin=15,
  chains=2,
  control=list(adapt_delta=0.99),
  save_warmup = FALSE
)

saveRDS(fit, file='results/STAN_model_test_Bspline_factorcov_simdata_full.RDS')
