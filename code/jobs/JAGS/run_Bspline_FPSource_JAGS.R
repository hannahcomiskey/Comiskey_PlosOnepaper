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

# Splines ------------------------------------------------------------
B <- splines::bs(all_years, degree=3, knots=c(seq(min(all_years), max(all_years), by=4)),  intercept = TRUE) # creating the B-splines
K <- dim(B)[2]
H = K-1

## The required data ------------------------------
inputdata <- list(Y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                  X = all_years,
                  num_years = n_all_years,
                  num_data = nrow(logit.data),
                  num_knots = length(B$knots.k),
                  P_count = length(n_subnat),
                  M_count = length(n_method),
                  C_count = length(n_country),
                  B = B,
                  matchsubnat = match_subnat,
                  matchmethod = match_method,
                  matchcountry = match_country,
                  matchyears = match_years
                  )

## Parameters to look at ------------------------------
pars <- c("a",
          "a0",
          "beta_c",
          "sigma_a",
          "sigma_a0",
          "P",
          "xi_sq",
          "kappa_sq",
          "a_xi",
          "c_xi",
          "kappa_sq_B",
          "sigma_beta")


# ## Run the model ---------------------
mod <- jags(data=inputdata,
            parameters.to.save=pars,
            model.file = "model/JAGS/Bspline_JAGS_TG_SEmodel.txt",
            n.iter = 80000,  # total number of iterations per chain
            n.burnin = 10000,
            n.thin=35)

saveRDS(mod, file='results/JAGS/JAGS_mod_Bspline_TG.RDS')


