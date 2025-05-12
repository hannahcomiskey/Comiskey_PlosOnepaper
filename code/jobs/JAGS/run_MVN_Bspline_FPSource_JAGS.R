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
inputdata <- list(Y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                  X = all_years,
                  num_years = n_all_years,
                  num_data = nrow(logit.data),
                  num_knots = K,
                  P_count = length(n_subnat),
                  M_count = length(n_method),
                  C_count = length(n_country),
                  B = B.ik,
                  Omega = (M+1)*diag(M),
                  matchsubnat = match_subnat,
                  matchmethod = match_method,
                  matchcountry = match_country,
                  matchyears = match_years
                  )
## Parameters to look at ------------------------------
pars <- c("a",
          "beta_c",
          "a0",
          "err_a0",
          "inv.Sigma.a0",
          "err_a",
          "inv.Sigma.a",
          "err_betac",
          "inv.Sigma.betac",
          "P",
          "Y_sim")


# ## Run the model ---------------------
mod <- jags(data=inputdata,
            parameters.to.save=pars,
            model.file = "model/JAGS/Bspline_MVN_JAGS_SEmodel.txt",
            n.iter = 80000,  # total number of iterations per chain
            n.burnin = 10000,
            n.thin=35)

saveRDS(mod, file='results/JAGS/JAGS_model_4country_MVN_Bspline.RDS')


