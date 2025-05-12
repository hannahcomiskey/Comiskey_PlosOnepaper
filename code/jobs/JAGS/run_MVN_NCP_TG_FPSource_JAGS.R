library(rjags)
library(R2jags)
library(tidyverse)
library(tidybayes)
library(parallel)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")


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
                  Omega = (M+1)*diag(M),
                  matchcountry= index_country_subnat_tbl$index_country,
                  matchsubnat = FP_source_data_wide$index_subnat,
                  matchmethod = FP_source_data_wide$index_method,
                  matchyears = FP_source_data_wide$index_year)

## Parameters to look at ------------------------------
pars <- c("alpha_pms", # required for P
          "delta.k",
          "alpha_cms",
          "beta.k",
          "err_alpha_cms",
          "err_alpha_pms",
          "var_delta",
          "P",
          "xi_sq",
          "kappa_sq",
          "a_xi",
          "c_xi",
          "kappa_sq_B")


# ## Run the model ---------------------
mod <- jags(data=inputdata,
            parameters.to.save=pars,
            model.file = "model/JAGS/ncp_MVN_TG_SEmodel_JAGS.txt",
            n.iter = 80000,         # total number of iterations per chain
            n.burnin = 10000,
            n.thin=35)

saveRDS(mod, file='results/JAGS/JAGS_model_MVN_NCP_TG_kstar.RDS')


