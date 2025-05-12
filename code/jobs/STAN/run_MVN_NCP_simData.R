library(rstan)
library(tidyverse)
library(tidybayes)
source('code/load_functions.R')

# Source simulated data --------------------------------------
options(mc.cores = parallel::detectCores())

load("data/simulated_data/simulated_data_MVN_kstar_Kenya_new.RData")

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
all_years <- seq(1:max(P_sim_df_sample$index_year))
B <- bs_bbase_precise(all_years)
Bik <- B$B.ik
K <-dim(Bik)[2]
H = K-1
kstar = B$Kstar


# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- rep(1, 6)

n_all_years <- length(all_years)
M_count = 5


Bik_array <- array(NA, dim=c(P,n_all_years,K))
for(i in 1:P) {
  Bik_array[i,,] <- Bik
  
}

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
                  matchyears = simmatchyears,
                  scale_global=scale_global)

## Parameters to look at ------------------------------
pars <- c("alpha_pms", # required for P
          "delta_k",
          "beta_c",
          "beta_k",
          "sigmabeta_Omega",
          "sigmabeta_tau",
          "sigmaalpha_Omega",
          "sigmaalpha_tau",
          "sigma_delta",
          "P")

# Run stan model ------------------

fit <- stan(
  file = "model/STAN/kstar_spline_public_private_NCP_logit_MVN_alpha_N_delta.stan",  # Stan program
  data = inputdata,    # named list of data
  pars = pars,
  iter = 10000,         # total number of iterations per chain
  warmup = 2000,
  thin=4,
  chains=3,
  save_warmup = FALSE,
  control=list(adapt_delta=0.99)
)


saveRDS(fit, file='results/STAN_model_kstar_NCP_MVN_SimKenya.RDS')
