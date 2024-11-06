library(rstan)
library(tidyverse)
library(tidybayes)

# Source code --------------------------------------
# source("read_data.R")
# source("prep_data.R")

# Read in SE data
subnat_FPsource_data <- readRDS("subnat_bivar_SE_source_data_20.RDS") # All countries

subnat_FPsource_data <- subnat_FPsource_data %>%
  ungroup() %>%
  dplyr::mutate(Region = stringr::str_to_title(Region))

# Updated cleaning approach -------------------
subnat_FPsource_data <- subnat_FPsource_data %>%
  dplyr::filter(Region!="NA")

# Filter where the sample size is smaller than 2 people across all three sectors ---------

FP_source_data_wide <- subnat_FPsource_data %>% # Proportion data
  dplyr::ungroup() %>%
  dplyr::rename(Public.SE = se.Public, Private.SE = se.Private) %>%
  dplyr::select(Country, Region, Method,  average_year, Public, Private, Public.SE, Private.SE, Public_n, Private_n) %>%
  dplyr::arrange(Country) %>%
  dplyr::filter(Public_n>=10 | Private_n>=10)

# Make sure proportions add to 1 ----------
FP_source_data_wide <- FP_source_data_wide %>%
  dplyr::rowwise() %>%
  dplyr::mutate(check_total = sum(Public, Private, na.rm = TRUE))

col_index <- which(colnames(FP_source_data_wide)=="Public")-1 # column index before CM column, as CM=1
# When check_Total=1, replace missing values with 0 ----------
for (i in 1:nrow(FP_source_data_wide)) {
  if(FP_source_data_wide$check_total[i]>0.99) {
    na_cols <- which(is.na(FP_source_data_wide[i, c("Public", "Private")])==TRUE) # NA values
    FP_source_data_wide[i, na_cols+col_index] <- as.list(rep(0, length(na_cols)))
  }
}

# Transform exactly 1 and 0 values away from boundary using lemon-squeezer approach ---------
FP_source_data_wide <- FP_source_data_wide %>%
  dplyr::mutate(Public = (Public*(nrow(FP_source_data_wide)-1)+0.5)/nrow(FP_source_data_wide)) %>%
  dplyr::mutate(Private = (Private*(nrow(FP_source_data_wide)-1)+0.5)/nrow(FP_source_data_wide)) %>%
  dplyr::select(Country, Region, Method, average_year, Public, Private, Public.SE, Private.SE, Public_n, Private_n, check_total) #, count_NA, remainder)

# Clean SE values --------------
FP_source_data_wide$count_SE.NA <- rowSums(is.na(FP_source_data_wide %>% dplyr::select(Public.SE, Private.SE))) # count NAs

SE_source_data_wide_norm <- FP_source_data_wide %>% dplyr::filter(count_SE.NA==0 & Private.SE>=0.005 & Public.SE>=0.005) # Normal obs. No action needed.
SE_source_data_wide_X <- FP_source_data_wide %>% dplyr::filter(Private.SE<0.005 | Public.SE<0.005) # Get obs with two missing sectors

# Replacing SE=0 and no NAs: Leontine's suggestion
SE_source_data_wide_X <- SE_source_data_wide_X %>%
  dplyr::filter(Public_n>=20 | Private_n>=20) # Remove small sample sizes (DHS has 10 units sampled per cluster as min., 20 as average)
col_index <- which(colnames(SE_source_data_wide_X)=="Public.SE")-1 # column index before CM column, as CM=1
DEFT_data <- readxl::read_xlsx("DEFT_DHS_database.xlsx") %>%
  rename(average_year = Year)
SE_source_data_wide_X <- SE_source_data_wide_X %>% left_join(DEFT_data)

# https://onlinestatbook.com/2/sampling_distributions/samp_dist_p.html

for(i in 1:nrow(SE_source_data_wide_X)) {
  num.SE0 <- which(SE_source_data_wide_X[i,c("Public.SE", "Private.SE")]<0.005)
  num.SEna <- which(is.na(SE_source_data_wide_X[i,c("Public.SE", "Private.SE")])==TRUE)
  DEFT <- ifelse(is.na(SE_source_data_wide_X$DEFT[i])==TRUE, 1.5, SE_source_data_wide_X$DEFT[i])
  N1 <- sum(SE_source_data_wide_X[i, c('Public_n', 'Private_n')], na.rm=TRUE) # Number of women surveyed
  phat <- 0.5/(N1+1) # Posterior mean of p under Jefferys prior for true prevalence of 0s.
  SE.hat <- sqrt((phat*(1-phat))/N1)
  SE_source_data_wide_X[i,c(col_index+num.SEna,col_index+num.SE0)] <- SE.hat*DEFT
}

FP_source_data_wide <- bind_rows(SE_source_data_wide_norm, SE_source_data_wide_X) # Put data back together again

# Remove proportions with two sectors still missing
FP_source_data_wide <- FP_source_data_wide %>%
  filter(is.na(Public)==FALSE)

FP_source_data_wide <- FP_source_data_wide %>% dplyr::arrange(Country, Region, Method, average_year)

# Replace issues with Burkina Faso names
FP_source_tmp <- FP_source_data_wide %>%
  dplyr::filter(Country=="Burkina Faso")  %>%
  dplyr::mutate(Region = case_when(Region == "North" ~ "Nord",
                                   Region == "East" ~ "Est",
                                   Region == "West" ~ "Centre-Ouest",
                                   Region ==  "Central/South" ~ "Centre-Sud",
                                   Region == "Ouagadougou" ~ "Centre Including Ouagadougou",
                                   TRUE ~ as.character(Region)))

FP_source_data_wide <- FP_source_data_wide %>%
  dplyr::filter(Country!="Burkina Faso")

FP_source_data_wide <- FP_source_data_wide %>%
  merge(FP_source_tmp, all = TRUE)

country_subnat_tbl <- FP_source_data_wide %>%
  group_by(Country, Region) %>%
  dplyr::select(Country, Region) %>%
  distinct() # table of country and regions (repeats in region names)

country_subnat_tbl <- country_subnat_tbl %>%
  ungroup() %>%
  mutate(index_subnat = 1:nrow(country_subnat_tbl))

country_tbl <- FP_source_data_wide %>%
  group_by(Country) %>%
  dplyr::select(Country) %>%
  distinct() # table of country and regions (repeats in region names)

country_tbl <- country_tbl %>%
  ungroup() %>%
  mutate(index_country = 1:nrow(country_tbl))

n_method <- c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills" ) # As per the method correlation matrix

method_index_table <- tibble(Method = n_method, index_method = 1:length(n_method))

n_country <- unique(country_subnat_tbl$Country)
n_subnat <- country_subnat_tbl$Region

FP_source_data_wide <- left_join(FP_source_data_wide, country_subnat_tbl)
FP_source_data_wide <- left_join(FP_source_data_wide, country_tbl)
FP_source_data_wide <- left_join(FP_source_data_wide, method_index_table)

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
n_sector <- c("Public", "Private") # Names of categories
n_obs <- nrow(FP_source_data_wide) # Total number of observations
year_seq <- seq(min(t_seq_2),max(t_seq_2), by=1)
n_years <- length(year_seq)

# Find the observation year indexes in the prediction years
country_index_tbl <- FP_source_data_wide %>%
  group_by(Country, index_country) %>%
  dplyr::select(Country, index_country) %>%
  distinct() # table of country and regions (repeats in region names)

index_country_subnat_tbl <- FP_source_data_wide %>%
  group_by(index_country, index_subnat) %>%
  dplyr::select(Region, index_subnat, Country, index_country) %>% #Super_region, index_superregion) %>%
  distinct() # table of country and regions (repeats in region names)

count_provinces <- index_country_subnat_tbl %>% # count number of provinces in each country
  group_by(index_country) %>%
  count(index_country) %>%
  dplyr::select(index_country, n) %>%
  distinct()

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

# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)

# Get logit of parameters and variance --------------------
mydata <- FP_source_data_wide[,c("Public", "Public.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))

# # testing splines ------------------------------------------------------------
B <- splines::bs(all_years, df=10, degree=3, intercept = FALSE)
K <-dim(B)[2]
B.ik <- B
D.hk <- diff(diag(K), diff = 1) # first order difference matrix (h = k-1)
Q.kh <- t(D.hk)%*%solve(D.hk%*%t(D.hk))
Zih <- B.ik%*%Q.kh
H <- dim(Zih)[2]

# Set up model inputs ----------------------------------------------------------
D=2
M_count = length(n_method)
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
pars <- c("L_t",
          "L_d",
          "psi",
          "sigma_alpha",
          "sigma_beta",
          "delta_k",
          "beta_c",
          "alpha_pms")

# Run stan model ---------------------------------------------------------------

fit <- stan(
  file = "zih_deltak_public_private_logit_SE_2.stan",  # Stan program
  data = inputdata,    # named list of data
  pars = pars,
  iter = 80000,         # total number of iterations per chain
  warmup = 10000,
  thin=35,
  chains=3,
  control=list(adapt_delta=0.99),
  save_warmup = FALSE
)

saveRDS(fit, file='STAN_model_zih_deltak_factorcov_Q_obsdata.RDS')
