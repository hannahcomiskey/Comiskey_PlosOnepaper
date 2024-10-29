library(rstan)
library(tidyverse)
library(tidybayes)
library(ggcorrplot)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/read_in_subnational_SEdata.R")

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


res <- cor(beta_df[,-1], use = "pairwise")
round(res, 2)

ggcorrplot(res, hc.order = TRUE, title= "Most recently observed public sector levels at the country level", type = "lower", lab = TRUE)
ggsave(filename = "estimated_correlations_countrylevel_intercepts.pdf", path = "visualisations/", height=12, width=15)

# Province-level intercepts
alpha_df <- FP_source_data_wide %>% 
  ungroup() %>%
  group_by(Country, Region, Method) %>%
  filter(average_year == max(average_year)) %>%
  ungroup() %>%
  select(Country, Region, Method, Public) %>%
  mutate(logit.Public = log(Public/(1-Public))) %>%
  select(!Public) %>%
  pivot_wider(names_from = Method, values_from = logit.Public)


res_alpha <- cor(alpha_df[,-c(1,2)], use = "pairwise")
round(res_alpha, 2)

ggcorrplot(res_alpha, hc.order = TRUE, title= "Most recently observed public sector levels at the provincial level", type = "lower", lab = TRUE)
ggsave(filename = "estimated_correlations_provincelevel_intercepts.pdf", path = "visualisations/", height=12, width=15)


# Investigating the constant variance assumption for the delta parameters

# Get the deviations away from alpha for each province, method
n_method <- c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills" ) # As per the method correlation matrix

logit_df <- FP_source_data_wide %>% 
  ungroup() %>%
  select(Country, Region, Method, average_year, Public) %>%
  mutate(logit.Public = log(Public/(1-Public))) %>%
  select(!Public) %>%
  pivot_wider(names_from = Method, values_from = logit.Public)

AD_df <- tibble(Country = NA, Region=NA , average_year=NA,  `Female Sterilization`=NA, Implants=NA, Injectables=NA, IUD=NA, `OC Pills`=NA) 

for(i in 1:nrow(logit_df)) {
  alpha_intercept <- alpha_df %>%
    filter(Country==logit_df$Country[i] & Region==logit_df$Region[i]) %>%
    select(all_of(n_method)) %>%
    unlist() %>%
    as.vector()
  
  tmp_logit_df <- logit_df[i, n_method] %>% unlist() %>% as.vector
  
  diff <- tmp_logit_df-alpha_intercept
  
  AD_df <- AD_df %>%
    add_row(Country=logit_df$Country[i], 
            Region=logit_df$Region[i],
            average_year=logit_df$average_year[i],
            `Female Sterilization`=diff[1], 
            Implants=diff[2], 
            Injectables=diff[3], 
            IUD=diff[4], 
            `OC Pills`=diff[5])
    
}

AD_df <- AD_df %>%
  filter(is.na(Country)==FALSE) 

AD_df %>%
  pivot_longer(cols=c( `Female Sterilization`, Implants, Injectables, IUD, `OC Pills`), names_to = 'Method', values_to = 'Difference') %>%
  ggplot() +
  geom_point(aes(x=average_year, y=Difference, colour=Method)) +
  geom_hline(yintercept = 0, lty='dashed')+
  theme(legend.position = 'bottom') 
ggsave(filename = "difference_most_recent_level.pdf", path = "visualisations/", height=12, width=15)

res_AD <- cor(AD_df[,-c(1,2,3)], use = "pairwise")
round(res_AD, 2)

ggcorrplot(res_AD, hc.order = TRUE, title= "Differences from most recently observed level at the province level", type = "lower", lab = TRUE)
ggsave(filename = "estimated_correlations_provincelevel_differences.pdf", path = "visualisations/", height=12, width=15)

# Correlations of rates of change
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

res_delta <- cor(delta_df[,c(9:13)], use = "pairwise")
round(res_delta, 2)

ggcorrplot(res_delta, hc.order = TRUE, title= "Rates of change in observed levels at the province level", type = "lower", lab = TRUE)
ggsave(filename = "estimated_correlations_provincelevel_ratesofchange.pdf", path = "visualisations/", height=12, width=15)


delta_df %>%
  pivot_longer(cols=c(lag.FS, lag.Implants, lag.Injectables, lag.IUD, lag.OCP), names_to = 'Method', values_to = 'delta') %>%
  ggplot() +
  geom_point(aes(x=average_year, y=delta, colour=Method)) +
  geom_hline(yintercept = 0, lty='dashed')+
  theme(legend.position = 'bottom') 
ggsave(filename = "delta_rates_of_change_over_time.pdf", path = "visualisations/", height=12, width=15)
