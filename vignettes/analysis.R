## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE
)

## ----subnational data---------------------------------------------------------
library(dplyr)
library(rjags)
library(R2jags)
library(tidyr)
library(ggplot2)
library(mcmsector)

glimpse(subnat_bivar_data)

## ----prepare data, echo=FALSE-------------------------------------------------

data <- clean_fp_source_data(mcmsector::subnat_bivar_data, deft_lookup = mcmsector::DEFT_DHS_database)

inputs <- run_preprocessing_pipeline(raw_df = data, 
                                     area_classification = mcmsector::Country_and_area_classification_inclFP2020,
                                     methods = c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills"),
                                     all_years = seq(1990, 2030.5, by = 0.5),
                                     nseg = 10)



## ----run JAGS model, echo=FALSE-----------------------------------------------
mydata <- inputs$data[,c("Public", "Public.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))
n_method = c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills")
M =  length(n_method)


## -----------------------------------------------------------------------------
## The required data ------------------------------
inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                  Bik = inputs$splines$B_ik,
                  n_years = inputs$meta$n_years,
                  n_obs = nrow(logit.data),
                  K = inputs$splines$K,
                  H = inputs$splines$H,
                  kstar = inputs$splines$Kstar,
                  C_count = length(inputs$meta$n_country),
                  P_count = length(inputs$meta$n_subnat),
                  M_count = length(c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills")),
                  Omega = diag(M) * 0.01 + diag(M),
                  matchcountry= inputs$meta$index_country_subnat_tbl$index_country,
                  matchsubnat = inputs$data$index_subnat,
                  matchmethod = inputs$data$index_method,
                  matchyears = inputs$data$index_year)

