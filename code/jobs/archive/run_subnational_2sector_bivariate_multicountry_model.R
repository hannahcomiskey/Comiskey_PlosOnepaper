library(R2jags)
library(runjags)
library(tidyverse)
library(tidybayes)
set.seed(1209)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")

# start parallel runs, save results in steps
library(foreach)
library(doMC)
registerDoMC()

n_iter = 80000
n_burnin = 10000
n_thin = 35
n_chains = 2

main_path = "results/global_subnat/bivar_2sector/"

# Get global correlations ------------------------------
estimated_rho_matrix <- readRDS("data/estimated_global_subnational_correlations.RDS")
estimated_rho_matrix <- estimated_rho_matrix %>%
  select(row, column, public_cor)
my_SE_rho_matrix <- estimated_rho_matrix %>%
  dplyr::select(public_cor) %>% 
  unlist() %>%
  as.vector()

# Get logit of parameters and variance --------------------
logit.data <- FP_source_data_wide[,c("Public", "Private")] %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Private = log(Private/(1-Private)))

logit.varcov_array <- readRDS("data/SE_source_data_VARCOV_updatedforSE_logit_bivar_2022.RDS") # includes RW AND SL


## The required JAGS data ------------------------------
jags_data <- list(y = logit.data[,c("logit.Public", "logit.Private")],
                  Sigma_y = logit.varcov_array,
                  rho = my_SE_rho_matrix, # correlations
                  kstar = Kstar,
                  B.ik = B.ik,
                  n_years = n_all_years,
                  n_obs = n_obs,
                  K = K,
                  H = H,
                  C_count = length(n_country),
                  P_count = length(n_subnat),
                  M_count = length(n_method),
                  matchsubnat = match_subnat,
                  matchcountry = match_country,
                  matchmethod = match_method,
                  matchyears = match_years
)

## Parameters to look at ------------------------------
jags_pars <- c("P",
              "alpha_pms", # required for P
              "alpha_cms", # required for single-country model
              "sigma_alpha_pms", 
              "sigma_alpha_cms",
              "sd_delta") # required for single-country model
              # "inv.sigma_delta",) # required for single-country model
#"delta.k") # required for correlations


## Do chains separately ------------------------------
foreach(chain=1:n_chains) %dopar% {
  set.seed(chain*1239)
  mod <- R2jags::jags(data = jags_data,
                      parameters.to.save = jags_pars,
                      model.file = "model/globalrun_Bspline_bivariate_logitNormal_2sector_model.txt",
                      n.chains = 1,
                      n.burnin = n_burnin,
                      n.iter = n_iter,
                      n.thin = n_thin,
                      working.directory = getwd()
  )
  ########### REMOVING SOME DATA FROM FIT OBJECT TO REDUCE FILE SIZE ####################################
  
  mod[names(mod) %in% c("model")] <- NA
  mod$BUGSoutput[names(mod$BUGSoutput) %in% c("sims.array", "sims.matrix", "mean", "sd")] <- NA
  saveRDS(mod, paste0(main_path,chain, "chain.rds"))
  print(paste("MCMC results for chain ", chain, "complete"))
} # end chains

gc()
chain=1
mod <- readRDS(paste0(main_path,chain, "chain.rds"))
for (chain in 2:n_chains) {
  mod_for_one_chain <- readRDS(paste0(main_path,chain, "chain.rds"))
  mod$BUGSoutput$sims.array <- mod$BUGSoutput$sims.array %>% abind::abind(mod_for_one_chain$BUGSoutput$sims.array, along = 2)
}
# now we need to hack the fit object such that it has correct meta data (as the original git object had meta data for just 1 chain)
mod$BUGSoutput$n.chains <- n_chains

## Save model ------------------------------
saveRDS(mod, paste0(main_path, "subnational_model_2sector_bivariate.RDS"))
