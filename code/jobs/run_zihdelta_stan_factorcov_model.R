library(rstan)
library(tidyverse)
library(tidybayes)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/read_in_subnational_SEdata.R")

FP_source_data_wide <- FP_source_data_wide %>% filter(Country=="Kenya")

area_classification <- mcmsupply::Country_and_area_classification %>% 
  dplyr::select(`Country or area`, Region) %>% 
  dplyr::rename(Country = `Country or area`) %>%
  dplyr::rename(Super_region = Region) %>%
  dplyr::mutate(Country = case_when(Country== "Bolivia (Plurinational State of)" ~ "Bolivia",
                                    Country== "Republic of Moldova" ~ "Moldova",
                                    Country== "Viet Nam" ~ "Vietnam",
                                    Country== "Kyrgyzstan" ~ "Kyrgyz Republic",
                                    .default = as.character(Country)))

FP_source_data_wide <- FP_source_data_wide %>%
  left_join(area_classification) 

#################################################
# Adding indexes to  data ----------------------
#################################################

country_subnat_tbl <- FP_source_data_wide %>% 
  group_by(Country, Region) %>% 
  dplyr::select(Country, Region) %>% 
  distinct() # table of country and regions (repeats in region names)
n_method <- c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills" ) # As per the method correlation matrix
n_country <- unique(country_subnat_tbl$Country)
n_subnat <- country_subnat_tbl$Region

FP_source_data_wide <- subnat_index_fun(FP_source_data_wide, n_subnat, country_subnat_tbl$Country)
FP_source_data_wide <- country_index_fun(FP_source_data_wide, n_country)
FP_source_data_wide <- method_index_fun(FP_source_data_wide, n_method)
# FP_source_data_wide <- superregion_index_fun(FP_source_data_wide, unique(FP_source_data_wide$Super_region))

# Time indexing - important for splines -----------------------------------
all_years <- seq(from = 1990, to = 2023.5, by=0.5) # shorter due to memory issues
n_all_years <- length(all_years)

FP_source_data_wide <- FP_source_data_wide %>%
  mutate(index_year = match(average_year,all_years)) 

#################################################
# setup for JAGS data ---------------------------
#################################################

t_seq_2 <- floor(FP_source_data_wide$index_year) # Time sequence for countries
country_seq <- FP_source_data_wide$Country
n_regions <- unique(FP_source_data_wide$Super_region)
n_sector <- c("Public", "Commercial_medical", "Other") # Names of categories
n_obs <- nrow(FP_source_data_wide) # Total number of observations
year_seq <- seq(min(t_seq_2),max(t_seq_2), by=1)
n_years <- length(year_seq)

# Find the observation year indexes in the prediction years
country_index_tbl <- FP_source_data_wide %>% 
  group_by(Country, index_country) %>% 
  dplyr::select(Country, index_country) %>%
  distinct() # table of country and regions (repeats in region names)

# index_superregion_tbl <- FP_source_data_wide %>% 
#   group_by(index_superregion, index_country) %>% 
#   dplyr::select(Country, index_country, Super_region, index_superregion) %>%
#   distinct() 

index_country_subnat_tbl <- FP_source_data_wide %>% 
  group_by(index_country, index_subnat) %>% 
  dplyr::select(Region, index_subnat, Country, index_country) %>% #Super_region, index_superregion) %>%
  distinct() # table of country and regions (repeats in region names)

count_provinces <- index_country_subnat_tbl %>% # count number of provinces in each country
  group_by(index_country) %>% 
  count(index_country) %>%
  dplyr::select(index_country, n) %>%
  distinct()

# index_superregion_province_tbl <- FP_source_data_wide %>% 
#   group_by(index_country, index_subnat) %>% 
#   dplyr::select(Region, index_subnat, Country, index_country) %>% #Super_region, index_superregion) %>%
#   distinct() 

# match_superregion <- index_superregion_tbl$index_superregion

match_country <- index_country_subnat_tbl$index_country

match_years <- FP_source_data_wide$index_year

match_method <- FP_source_data_wide$index_method

match_subnat <-  FP_source_data_wide$index_subnat

# Get T_star and match_Tstar -----------------------------
T_star <- FP_source_data_wide %>%
  group_by(Country, Region) %>%
  dplyr::filter(index_year==max(index_year)) %>%
  dplyr::select(Country, Region, index_country, index_subnat, average_year, index_year) %>%
  arrange(index_subnat) %>%
  ungroup() %>%
  dplyr::select(index_country, index_subnat, average_year, index_year) %>%
  distinct()

#################################################
# Splines --------------------------------------
#################################################
B <- splines::bs(all_years, df=10, degree=3, intercept = FALSE)
K <-dim(B)[2]
B.ik <- B
D.hk <- diff(diag(K), diff = 1) # first order difference matrix (h = k-1)
Q.kh <- t(D.hk)%*%solve(D.hk%*%t(D.hk))
Zih <- B.ik%*%Q.kh 
H <- dim(Zih)[2]

options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 1, auto_write = TRUE)

# Get logit of parameters and variance --------------------
mydata <- FP_source_data_wide[,c("Public", "Public.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))

# Set up model inputs ----------------------------------------------------------
D=2
M_count = 5
OD_count  = D*(M_count-D)+ D*(D-1)/2

## The required data -----------------------------------------------------------
inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                  Zih = Zih,
                  intercept = rep(1, n_all_years),
                  n_years = n_all_years,
                  delta_mu = rep(0, 5),
                  beta_mu = rep(0, 5),
                  n_obs = nrow(logit.data),
                  K = K,
                  H = H,
                  D=2,
                  OD_count = OD_count,
                  C_count = length(n_country),
                  P_count = length(n_subnat),
                  M_count = length(n_method),
                  S_count = 2,
                  matchsubnat = match_subnat,
                  matchcountry = match_country,
                  matchmethod = match_method,
                  matchyears = match_years)


## Parameters to look at -------------------------------------------------------
pars <- c("L_Sigma",
          "Q",
          "beta_c",
          "alpha_pms",
          "sigma_alpha",
          "delta_k")

# Run stan model ---------------------------------------------------------------

fit <- stan(
  file = "model/STAN/zih_deltak_factorcov_public_private_logit.stan",  # Stan program
  data = inputdata,    # named list of data
  pars = pars,
  iter = 40000,         # total number of iterations per chain
  warmup = 10000,
  thin=15,
  chains=2,
  control=list(adapt_delta=0.99),
  save_warmup = FALSE
)

saveRDS(fit, file='results/STAN_model_zih_deltak_factorcov_simdata_Kenya.RDS')
