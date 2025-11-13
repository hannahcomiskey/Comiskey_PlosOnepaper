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

Tstar_df <- FP_source_data_wide %>% select(Country, Region, index_year) %>%
  group_by(Country, Region) %>%
  filter(index_year == max(index_year)) %>%
  distinct()

## The required data ------------------------------
inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                  n_years = n_all_years,
                  n_obs = nrow(logit.data),
                  C_count = length(n_country),
                  P_count = length(n_subnat),
                  M_count = length(n_method),
                  matchcountry= index_country_subnat_tbl$index_country,
                  matchsubnat = FP_source_data_wide$index_subnat,
                  matchmethod = FP_source_data_wide$index_method,
                  matchyears = FP_source_data_wide$index_year)

## Parameters to look at ------------------------------
pars <- c("sd_beta",
          "P")


# ## Run the model ---------------------
mod <- jags(data=inputdata,
            parameters.to.save=pars,
            model.file = "model/JAGS/ncp_MVN_SEmodel_RW.txt",
            n.iter = 120000,         # total number of iterations per chain
            n.burnin = 20000,
            n.thin=50)

saveRDS(mod, file='results/JAGS/JAGS_model_MVN_RW.RDS')


