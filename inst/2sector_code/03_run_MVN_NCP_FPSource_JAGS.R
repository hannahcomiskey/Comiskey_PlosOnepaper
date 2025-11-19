library(rjags)
library(R2jags)
library(tidyverse)
library(tidybayes)
library(parallel)

# Source code --------------------------------------
source("R/load_functions.R")
source("R/2sector_code/01_read_in_subnational_2sector_data.R") # Read in data
source("R/2sector_code/02_set_up_2sector_bivar_globalrunjags.R") # Data cleaning and set up model inputs


# Get logit of parameters and variance --------------------
mydata <- FP_source_data_wide[,c("Public", "Public.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))
M =  length(n_method)

## The required data ------------------------------
inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                  Bik = B.ik,
                  n_years = n_all_years,
                  n_obs = nrow(logit.data),
                  K = K,
                  H = H,
                  kstar = Kstar,
                  C_count = length(n_country),
                  P_count = length(n_subnat),
                  M_count = length(n_method),
                  Omega = diag(M) * 0.01 + diag(M),
                  matchcountry= index_country_subnat_tbl$index_country,
                  matchsubnat = FP_source_data_wide$index_subnat,
                  matchmethod = FP_source_data_wide$index_method,
                  matchyears = FP_source_data_wide$index_year)

## Parameters to look at ------------------------------
pars <- c( "alpha_pms", # required for P
          "inv.Sigma.alpha_cms",
          "inv.Sigma.alpha_pms",
          "sigma_delta",
           "beta.k")


# ## Run the model ---------------------
mod <- jags.parallel(data=inputdata,
            parameters.to.save=pars,
            model.file = "inst/model/ncp_MVN_SEmodel_JAGS.txt",
            n.iter = 80000,         # total number of iterations per chain
            n.burnin = 10000,
            n.thin=35)

# saveRDS(mod, file='results/JAGS/JAGS_model_MVN_NCP_kstar_SE.RDS')


