library(rstan)
library(tidyverse)
library(tidybayes)

# Source simulated data --------------------------------------
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
simmatchcountry <- rep(1, 6)

n_all_years <- length(all_years)
M_count = 5


Bik_array <- array(NA, dim=c(P,n_all_years,K))
for(i in 1:P) {
  Bik_array[i,,] <- Bik
  
}

# Source simulated data --------------------------------------
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# let p0=3, D=H and n=n_obs 
# See 3.12 Vehtari paper.
scale_global = 3/((H-3)*sqrt(nrow(logit.data)))

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
                  scale_global=scale_global
                  )

## Parameters to look at ------------------------------
pars <- c("alpha_pms", # required for P
          "delta_k",
          "beta_c",
          "beta_k",
          "sigmabeta_Omega",
          "sigmabeta_tau",
          'sigma_alpha',
          "lstar_sq",
          "tau",
          "P")

# Run stan model ------------------

fit <- stan(
  file = "model/STAN/ncp_reg_MVN_model.stan",  # Stan program
  data = inputdata,    # named list of data
  pars = pars,
  iter = 10000,         # total number of iterations per chain
  warmup = 2000,
  thin=4,
  chains=3,
  save_warmup = FALSE,
  control=list(adapt_delta=0.99)
)


saveRDS(fit, file='results/STAN_model_kstar_reg_NCP_MVN_alpha_SimKenya.RDS')
