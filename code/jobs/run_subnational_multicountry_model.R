library(mcmsupply)
library(R2jags)
library(runjags)
library(tidyverse)
library(tidybayes)
set.seed(1209)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/read_in_subnational_SEdata.R")
source("code/set_up_globalrunjags.R")

# start parallel runs, save results in steps
library(foreach)
library(doMC)
registerDoMC()

n_iter = 80000
n_burnin = 10000
n_thin = 35
n_chains = 2

main_path = "results/global_subnat/test_corr/"

# Get global correlations ------------------------------
estimated_rho_matrix <- readRDS("data/estimated_global_subnational_correlations.RDS")
estimated_rho_matrix <- estimated_rho_matrix %>%
  select(row, column, public_cor, private_cor)
my_SE_rho_matrix <- estimated_rho_matrix %>%
  dplyr::select(public_cor, private_cor)

# Get logit of parameters and variance --------------------
mydata <- FP_source_data_wide[,c("Public", "Public.SE", "Commercial_medical", "Commercial_medical.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.CM = log(Commercial_medical/(1-Commercial_medical)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var),
         logit.CM.Var = ((1/(Commercial_medical*(1-Commercial_medical)))^2)*Commercial_medical.SE^2,
         logit.CM.SE = sqrt(logit.CM.Var))

# # Get number of complete observations per sector ----------------------
# logit.CM_ratio_num_obs <- FP_source_data_wide %>%
#   rowwise() %>%
#   mutate(CM_ratio = Commercial_medical/(1-Public),
#          logit.CM = log(Commercial_medical/(1-Commercial_medical)),
#          logit.CM_ratio = log(CM_ratio/(1-CM_ratio))) %>%
#   group_by(Country, Method, Region) %>%
#   filter(average_year==max(average_year)) %>% 
#   select(Method, logit.CM_ratio) %>%
#   pivot_wider(names_from=Method, values_from = logit.CM_ratio) %>%
#   #drop_na() %>%
#   nrow()
# 
# 
# logit.public_num_obs <- FP_source_data_wide %>%
#   rowwise() %>%
#   mutate(logit.public = log(Public/(1-Public))) %>%
#   group_by(Country, Method, Region) %>%
#   filter(average_year==max(average_year)) %>% 
#   select(Method, logit.public) %>%
#   pivot_wider(names_from=Method, values_from = logit.public) %>%
#   #drop_na() %>%
#   nrow()
# 
# n_complete <- c(logit.public_num_obs, logit.CM_ratio_num_obs)

## The required JAGS data ------------------------------
jags_data <- list(y = logit.data[,c("logit.Public", "logit.CM")], # using total proportions as collapsing over sectors
                  se_prop = logit.data[,c("logit.Public.SE", "logit.CM.SE")],
                  rho = my_SE_rho_matrix, # correlations
                  kstar = Kstar,
                  B.ik = B.ik,
                  n_years = n_all_years,
                  n_obs = n_obs,
                  K = K,
                  H = H,
                  C_count = length(n_country),
                  P_count = length(n_subnat),
                  R_count = length(n_regions),
                  M_count = length(n_method),
                  matchsubnat = match_subnat,
                  matchcountry = match_country,
                  matchregion = match_superregion,
                  matchmethod = match_method,
                  matchyears = match_years,
                  Omega = diag(length(n_method)),
                  J = (length(n_method)+1),
                 sd.Z = sqrt(0.2) # as per paper
)

## Parameters to look at ------------------------------
jags_pars <- c("alpha_pms", # required for P
               "alpha_cms", # required for single-country model
               "Sigma.alpha_cms", 
               "Sigma.alpha_pms", # required for single-country model
               "rho.alpha_cms",
               "sd.Z",
               "beta.k") #, # required for P
              # "inv.sigma_delta",) # required for single-country model
#"delta.k") # required for correlations


## Do chains separately ------------------------------
foreach(chain=1:n_chains) %dopar% {
  set.seed(chain*1239)
  mod <- R2jags::jags(data = jags_data,
                      parameters.to.save = jags_pars,
                      model.file = "model/globalrun_Bspline_logitNormal_estRho_model.txt",
                      n.chains = 1,
                      n.burnin = n_burnin,
                      n.iter = n_iter,
                      n.thin = n_thin,
                      working.directory = getwd()
  )
  ########### REMOVING SOME DATA FROM FIT OBJECT TO REDUCE FILE SIZE ####################################
  
  mod[names(mod) %in% c("model")] <- NA
  mod$BUGSoutput[names(mod$BUGSoutput) %in% c("sims.list",
                                              "summary",
                                              "mean",
                                              "sd")] <- NA
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
saveRDS(mod, paste0(main_path, "subnational_model_all_mvn_alpha_4country_sdZ_perpaper.RDS"))
