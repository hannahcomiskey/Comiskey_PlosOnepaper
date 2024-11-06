library(rstan)
library(tidyverse)
library(tidybayes)

# Source code --------------------------------------
# source("code/load_functions.R")
# source("code/2sector_code/read_in_subnational_2sector_data.R")
# source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")

# Source code --------------------------------------
source("code/load_functions.R")
source("code/read_in_subnational_SEdata.R")
source("code/set_up_globalrunjags.R")

options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 1, auto_write = TRUE)

# Get logit of parameters and variance --------------------
mydata <- FP_source_data_wide[,c("Public", "Public.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))

# logit.data <- FP_source_data_wide[,c("Public", "Private")] %>%
#   rowwise() %>%
#   mutate(logit.Public = log(Public/(1-Public)),
#          logit.Private = log(Private/(1-Private)))
# 
# logit.varcov_array <- readRDS("data/SE_source_data_VARCOV_updatedforSE_logit_bivar_2022.RDS") # includes RW AND SL
# 
# varcov_vec <- list()
# for(i in 1:nrow(FP_source_data_wide)) {
#   varcov_vec[[i]] <- logit.varcov_array[,,i]
# }
# 
# varcov_vec <- varcov_vec %>% as.vector()

# # testing splines ---------------------------------
B <- splines::bs(all_years, df=10, degree=3, intercept = FALSE)
K <-dim(B)[2]
B.ik <- B

## using a different implementation for the penalisation so that we can put a prior directly on the coefficient differences
## Eilers (1999) proposed Z = BD'(DD')^(-1) (where D is the differencing matrix)
## then instead of y = B*beta, we have y = alpha + Z*delta
## where delta are the first order differences (if delta is 0 then splines stay flat)
## this works much better for convergence

D.hk <- diff(diag(K), diff = 1) # first order difference matrix (h = k-1)
Q.kh <- t(D.hk)%*%solve(D.hk%*%t(D.hk))
Zih <- B.ik%*%Q.kh 
H <- dim(Zih)[2]


## The required data ------------------------------
inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                  Zih = Zih,
                  intercept = rep(1, n_all_years),
                  n_years = n_all_years,
                  splinestart = rep(0,length(n_subnat)),
                  delta_mu = rep(0, 5),
                  beta_mu = rep(0, 5),
                  n_obs = n_obs,
                  K = K,
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
pars <- c("alpha_pms", # required for P
          "delta_k",
          "P",
          "y_tilde",
          "L_Sigma")

# Run stan model ------------------

fit <- stan(
  file = "model/STAN/CopyOfzih_deltak_public_private_logit.stan",  # Stan program
  data = inputdata,    # named list of data
  pars = pars,
  iter = 80000,         # total number of iterations per chain
  warmup = 10000,
  thin=35,
  chains=3
)


saveRDS(fit, file='results/STAN_model_test_zih_cholesky.RDS')
