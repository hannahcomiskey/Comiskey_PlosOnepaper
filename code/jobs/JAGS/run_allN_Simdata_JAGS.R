library(rjags)
library(R2jags)
library(tidyverse)
library(tidybayes)
library(parallel)
source('code/load_functions.R')

# Source simulated data --------------------------------------
load("data/simulated_data/simulated_data_all_N_kstar_Kenya_10years.RData")

# Get logit of parameters and variance -----------------------------------------
mydata <- P_sim_df_sample[,c("Public")] %>%
  mutate(Public.SE = 0.1) %>%
  rowwise() %>%
  mutate(Public = ifelse(Public < 0.0001 , 0, Public))

logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)))

# # testing splines ------------------------------------------------------------
all_years <-  -5:55 # -5:20 #
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
          "delta.k",
          "alpha_cms",
          "beta.k",
          "sigma_alpha_pms",
          "sigma_alpha_cms",
          "sigma_delta",
          "P",
          "y.sim")


# ## Run the model ---------------------
mod <- jags(data=inputdata,
                     parameters.to.save=pars,
                     model.file = "model/ncp_allN_model_JAGS.txt",
                     n.iter = 60000,         # total number of iterations per chain
                     n.burnin = 10000,
                     n.thin=25)

saveRDS(mod, file='results/JAGS_model_KenyaSim_allN_NCP_kstar_sum0_10years.RDS')


