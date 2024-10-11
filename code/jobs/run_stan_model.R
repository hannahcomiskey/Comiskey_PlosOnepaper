# install.packages("devtools")
#devtools::install_github("r-lib/conflicted")

library(rstan)
library(tidyverse)
library(tidybayes)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/read_in_subnational_SEdata.R")
source("code/set_up_globalrunjags.R")

options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 1, auto_write = TRUE)

#example(stan_model, package = "rstan", run.dontrun = TRUE)

# # Get global correlations ------------------------------
# estimated_rho_matrix <- readRDS("data/estimated_global_subnational_correlations.RDS")
# estimated_rho_matrix <- estimated_rho_matrix %>%
#   select(row, column, public_cor)
# my_SE_rho_matrix <- estimated_rho_matrix %>%
#   dplyr::select(public_cor)

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

# logit.data <- FP_source_data_wide[,c("Public", "Public.SE")] %>%
#   rowwise() %>%
#   mutate(logit.Public = log(Public/(1-Public))) #,
#          #logit.Private = log(Private/(1-Private)))
# 
# logit.varcov_array <- readRDS("data/SE_source_data_VARCOV_updatedforSE_logit_bivar_2022.RDS") # includes RW AND SL
# 
# #varcov_vec <- list()
# varcov_vec <- vector()
# for(i in 1:nrow(FP_source_data_wide)) {
#   #varcov_vec[[i]] <- logit.varcov_array[,,i]
#   varcov_vec[i] <- sqrt(diag(logit.varcov_array[,,i])[1])
# }


# # testing splines ---------------------------------
B <- splines::bs(all_years, df=10, degree=3, intercept = FALSE)

# ## Plot basis
# par(lwd = 3, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3, mfrow = c(1,1))
# plot(all_years, B[,1],  type= "l", col = 1, lwd = 1,
#      ylim = c(0,1),
#      ylab ="Basis Function")
# for (k in 1:ncol(B)){
#   lines(all_years, B[,k], type= "l", col = k, lwd = 1)
# }

Bik_vec <- list()
for(i in 1:length(n_subnat)) {
  Bik_vec[[i]] <-B
}

Bik_vec <- Bik_vec %>% as.vector()
  
## The required data ------------------------------
inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                  #Sigma_Y = varcov_vec,
                  kstar = 10, #Kstar,
                  Bik = B,
                  intercept = rep(1, n_all_years),
                  n_years = n_all_years,
                  splinestart = rep(0,length(n_subnat)),
                  beta_mu = rep(0, 5),
                  delta_mu = rep(0, 5),
                  n_obs = n_obs,
                  K = ncol(B),
                  H = H,
                  C_count = length(n_country),
                  P_count = length(n_subnat),
                  M_count = length(n_method),
                  S_count = 2,
                  matchsubnat = match_subnat,
                  matchcountry = match_country,
                  matchmethod = match_method,
                  matchyears = match_years
)

## Parameters to look at ------------------------------
# pars <- c("alpha_pms", # required for P
#           "z")

# Run stan model ------------------

fit <- stan(
  file = "model/STAN/spline_model.stan",  # Stan program
  data = inputdata,    # named list of data
  iter = 200,         # total number of iterations per chain
  chains=1
  )

saveRDS(fit, file='results/STAN_model_test.RDS')
