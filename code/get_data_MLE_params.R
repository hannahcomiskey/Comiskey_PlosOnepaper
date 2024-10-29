library(tidyverse)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/2sector_code/read_in_subnational_2sector_data.R")

# Get intercepts for each province-method combo 
# Most recently observed levels
# Country-level intercepts
beta_df <- FP_source_data_wide %>% 
  ungroup() %>%
  group_by(Country, Method) %>%
  filter(average_year == max(average_year)) %>%
  ungroup() %>%
  select(Country, Region, Method, Public) %>%
  mutate(logit.Public = log(Public/(1-Public))) %>%
  select(!Public) %>%
  pivot_wider(names_from = Method, values_from = logit.Public) %>%
  select(!Region) %>%
  ungroup() %>%
  group_by(Country) %>%
  mutate(mean.OCP = mean(`OC Pills`, na.rm=TRUE),
         mean.Injectables = mean(Injectables, na.rm=TRUE),
         mean.IUD = mean(IUD, na.rm=TRUE),
         mean.FS = mean(`Female Sterilization`, na.rm=TRUE),
         mean.Implants = mean(Implants, na.rm=TRUE)) %>%
  select(Country, mean.FS, mean.Implants, mean.Injectables, mean.IUD, mean.OCP) %>%
  distinct()

saveRDS(beta_df, 'data/simulated_data/mle_mean_beta.RDS')

# Varcov
head(FP_source_data_wide)
n_method <- c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills" ) # As per the method correlation matrix

alpha_df <- FP_source_data_wide %>% 
  ungroup() %>%
  group_by(Country, Region, Method) %>%
  filter(average_year == max(average_year)) %>%
  ungroup() %>%
  select(Country, Region, Method, Public) %>%
  mutate(logit.Public = log(Public/(1-Public))) %>%
  select(!Public) %>%
  pivot_wider(names_from = Method, values_from = logit.Public)

alpha_df[which(complete.cases(alpha_df[,c(3:7)])==TRUE),]

cov_alpha <- as_tibble(cov(alpha_df[,c(3:7)], use='complete.obs')) %>%
  mutate(Method = rownames(cov(alpha_df[,c(3:7)], use='complete.obs'))) %>%
  select(Method, all_of(n_method)) %>%
  arrange(factor(Method, levels = n_method))

mle_cov_alpha <- as.matrix(cov_alpha[,2:6])
saveRDS(mle_cov_alpha, 'data/simulated_data/mle_cov_alpha.RDS')

# Get rates of change for each province-method combo 
# Rates of change from most recently observed level 

delta_df <- FP_source_data_wide %>% 
  ungroup() %>%
  group_by(Country, Region, Method) %>%
  arrange(Country, Region, Method, average_year) %>%
  ungroup() %>%
  select(Country, Region, Method, average_year, Public) %>%
  mutate(logit.Public = log(Public/(1-Public))) %>%
  select(!Public) %>%
  pivot_wider(names_from = Method, values_from = logit.Public) %>%
  group_by(Country, Region) %>%
  arrange(Country, Region,average_year) %>%
  mutate(lag.OCP = `OC Pills`-lag(`OC Pills`),
         lag.Injectables = Injectables-lag(Injectables),
         lag.IUD = IUD-lag(IUD),
         lag.FS = `Female Sterilization`-lag(`Female Sterilization`),
         lag.Implants = Implants-lag(Implants))

delta_df[which(complete.cases(delta_df[,c(9:13)])==TRUE),]

var.OCP <- var(delta_df$lag.OCP, na.rm=TRUE)
var.Inj <- var(delta_df$lag.Injectables, na.rm=TRUE)
var.IUD <- var(delta_df$lag.IUD, na.rm=TRUE)
var.FS <- var(delta_df$lag.FS, na.rm=TRUE)
var.Imp <- var(delta_df$lag.Implants, na.rm=TRUE)

n_method <- c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills" ) # As per the method correlation matrix

mle_var_delta <- c(var.FS, var.Imp, var.Inj, var.IUD, var.OCP)

saveRDS(mle_var_delta, 'data/simulated_data/mle_var_delta.RDS')



