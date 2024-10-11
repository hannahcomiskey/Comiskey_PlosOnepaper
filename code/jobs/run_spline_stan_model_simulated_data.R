library(rstan)
library(tidyverse)
library(tidybayes)

# Source simulated data --------------------------------------
options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 1, auto_write = TRUE)

load("data/simulated_data_new.RData")

# Get logit of parameters and variance --------------------
mydata <- P_sim_df_sample[,c("Public")] %>%
  mutate(Public.SE = 0.1) %>%
  rowwise() %>%
  mutate(Public = ifelse(Public < 0.0001 , 0, Public))

logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))

# # testing splines ---------------------------------
B <- splines::bs(all_years, df=10, degree=3, intercept = FALSE)
K <-dim(B)[2]
B.ik <- B

# Set up model inputs -----------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- matchcountry
n_all_years <- length(all_years)

## The required data ------------------------------
inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  # se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                  Bik = B.ik,
                  intercept = rep(1, n_all_years),
                  n_years = n_all_years,
                  beta_mu = rep(0, 5),
                  n_obs = nrow(logit.data),
                  K = K,
                  H = H,
                  C_count = C,
                  P_count = P,
                  M_count = M,
                  S_count = 2,
                  matchsubnat = simmatchsubnat,
                  matchcountry = simmatchcountry,
                  matchmethod = simmatchmethod,
                  matchyears = simmatchyears
)

## Parameters to look at ------------------------------
pars <- c("alpha_pms", # required for P
          "beta_c",
          "sigma_alpha",
          "sigma_beta",
          "a_k",
          "tau_a",
          "a_raw",
          "sigma_y",
          "P")

# Run stan model ------------------

fit <- stan(
  file = "model/STAN/penalised_Bspline_model.stan",  # Stan program
  data = inputdata,    # named list of data
  pars = pars,
  iter = 50000,         # total number of iterations per chain
  warmup = 10000,
  thin=20,
  chains=3,
  save_warmup = FALSE,
  control=list(adapt_delta=0.99)
)


saveRDS(fit, file='results/STAN_model_test_spline_simdata_sigmay.RDS')
