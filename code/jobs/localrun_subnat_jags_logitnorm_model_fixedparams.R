# Source code --------------------------------------
library(R2jags)
library(runjags)
library(tidyverse)
library(tidybayes)

# Load data ----------------------------------------
source("code/load_functions.R")
source("code/read_in_subnational_SEdata.R")
source("code/get_local_params.R")
source("code/set_up_jags_SE_data_refalpha.R")

# Name of country --------------------------------------
print(subset_C)

# start parallel runs, save results in steps
library(foreach)
library(doMC)
registerDoMC()

n_iter = 70000
n_burnin = 20000
n_thin = 25
n_chains = 2

# Create folders to save results in ---------------------------
tmp_res = "results/fixed_params/"
tmp_folders <- list.files(tmp_res)
if(subset_C %in% tmp_folders){
  main_path = paste0(tmp_res,subset_C,"/")
} else { 
  dir.create(paste0(tmp_res,subset_C,"/"))
  main_path = paste0(tmp_res,subset_C,"/")
}

# Get intercepts -------------------
myalpha_med <- alpha_cms_hat[,,subset_C] # Take out relevant region

# Get logit of parameters and variance --------------------
mydata <- FP_source_subset[,c("Public", "Public.SE", "Commercial_medical", "Commercial_medical.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.CM = log(Commercial_medical/(1-Commercial_medical)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var),
         logit.CM.Var = ((1/(Commercial_medical*(1-Commercial_medical)))^2)*Commercial_medical.SE^2,
         logit.CM.SE = sqrt(logit.CM.Var))

## The required JAGS data ------------------------------
jags_data <- list(y = logit.data[,c("logit.Public", "logit.CM")], # using total proportions as collapsing over sectors
                  se_prop = logit.data[,c("logit.Public.SE", "logit.CM.SE")],
                  alpha_cms_hat = myalpha_med,
                  tau_alpha_snms_hat = tau_alpha_pms_hat,
                  inv.sigma_delta = inv.sigma_delta_hat, # dwish on inverse
                  tstar = T_star$index_year,
                  kstar = Kstar,
                  B.ik = B.ik,
                  n_years = n_all_years,
                  n_obs = n_obs,
                  K = K,
                  H = H,
                  P_count = length(n_subnat),
                  M_count = length(n_method),
                  matchsubnat = match_subnat,
                  matchcountry = match_country,
                  matchmethod = match_method,
                  matchyears = match_years)

## Parameters to look at ------------------------------
jags_pars <- c("P",
               "alpha_snms",
               "delta.k")

## Do chains separately ------------------------------
foreach(chain=1:n_chains) %dopar% {
  set.seed(chain*1239)
  mod <- R2jags::jags(data = jags_data,
                      parameters.to.save = jags_pars,
                      model.file = "model/singlecountry_logitnorm_model_fixed_params.txt",
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
saveRDS(mod, paste0(main_path,"mod_fixedparms_",subset_C,".RDS"))

# Get local point summary ------------------
source("code/get_local_point_estimates.R")
