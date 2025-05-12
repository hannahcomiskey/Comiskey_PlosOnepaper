library(rjags)
library(R2jags)
library(tidyverse)
library(tidybayes)
library(parallel)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")

# Read in results from original model 
original_P_samps <- readRDS('results/global_subnat/bivar_2sector/JAGS_original_Psamps_df.RDS') %>%
  mutate(Source = 'Original')

# Get logit of parameters and variance --------------------
mydata <- FP_source_data_wide[,c("Public", "Public.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))
M =  length(n_method)

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))

# Read in JAGS results ---------------------------------------------------------
jags_N_p_samps <- readRDS('results/JAGS/P_samps/N_delta/JAGS_N_delta_Psamps_df.RDS') %>%
  mutate(Source = 'JAGS N delta model')

jags_MVN_p_samps <- readRDS('results/JAGS/P_samps/MVN_delta_chains/JAGS_MVN_delta_Psamps_df.RDS') %>%
  mutate(Source = 'JAGS MVN delta model')

jags_Bspline_N_p_samps <- readRDS('results/JAGS/P_samps/Bspline/JAGS_N_bspline_Psamps_df.RDS') %>%
  mutate(Source = 'JAGS N Bspline model')

jags_Bspline_MVN_p_samps <- readRDS('results/JAGS/P_samps/MVN_bspline/JAGS_MVN_bspline_Psamps_df.RDS') %>%
  mutate(Source = 'JAGS MVN Bspline model')

# Read in STAN results ---------------------------------------------------------
stan_N_p_samps <- readRDS('results/STAN/P_samps/STAN_model_4country_allN_NCP_Psamps_df.RDS') %>%
  mutate(Source = 'STAN N delta model')

stan_Bspline_N_p_samps <- readRDS('results/STAN/P_samps/STAN_model_4country_allN_NCP_Psamps_df.RDS') %>%
  mutate(Source = 'STAN N Bspline model')

# Get observed data 
P_df<- FP_source_data_wide %>%
  select(Country, Region, Method, index_year, average_year, Public, Private) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

P_SEdf<- FP_source_data_wide %>%
  select(Country, Region, Method, index_year, average_year, Public.SE, Private.SE) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public.SE, Private.SE), names_to = 'Sector', values_to = 'SE') %>%
  mutate(Sector = str_replace(Sector, '.SE', ''))

P_df <- left_join(P_df, P_SEdf) %>%
  rowwise() %>%
  mutate(lower_95 = Observed - 2*SE,
         upper_95 = Observed + 2*SE)

for(c in n_country) {
  ggplot() +
    geom_point(data = P_df %>% filter(Country==c ), aes(x=average_year, y=Observed, colour=Sector, pch=Sector)) +
    geom_errorbar(data = P_df %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax=upper_95, colour=Sector)) +
    geom_line(data = jags_N_p_samps %>% filter(Country==c ), aes(x=average_year, y=Mean, colour=Sector, lty=Source)) +
    geom_ribbon(data = jags_N_p_samps %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
    geom_line(data = stan_N_p_samps %>% filter(Country==c ) , aes(x=average_year, y=Mean, colour=Sector, lty=Source)) +
    geom_ribbon(data = stan_N_p_samps %>% filter(Country==c ), aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
    theme(title = element_text(size=10), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
    facet_wrap(~interaction(Method, Region), ncol = 5)
  ggsave(filename = paste0('compare_p_plot_',c,'_N_kstar_JAGS_STAN.pdf'), path = 'visualisations/N_delta/', height=12, width=15) 
  
}

n_sector = c('Public', 'Private')

# Plot means vs observed values
for(c in n_country) {
  for(s in n_sector) {
    print(c)
    print(c)
    ggplot() +
      geom_point(data = P_df %>% filter(Country==c & Sector==s), aes(x=average_year, y=Observed, pch=Sector)) +
      geom_errorbar(data = P_df %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax=upper_95)) +
      geom_line(data = original_P_samps %>% filter(Country==c & Sector==s), aes(x=average_year, y=Mean, colour=Source, lty=Source)) +
      geom_ribbon(data = original_P_samps %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax = upper_95,  lty=Source, fill=Source), alpha=0.2) +
      geom_line(data = jags_N_p_samps %>% filter(Country==c & Sector==s), aes(x=average_year, y=Mean, colour=Source, lty=Source)) +
      geom_ribbon(data = jags_N_p_samps %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax = upper_95, lty=Source, fill=Source), alpha=0.2) +
      # geom_line(data = stan_N_p_samps %>% filter(Country==c & Sector==s) , aes(x=average_year, y=Mean, colour=Source, lty=Source)) +
      # geom_ribbon(data = stan_N_p_samps %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax = upper_95, lty=Source, fill=Source), alpha=0.2) +
      theme(title = element_text(size=10), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
      facet_wrap(~interaction(Method, Region), ncol = 5)
    ggsave(filename = paste0('compare_p_plot_',c,'_',s,'_N_kstar_JAGS_original.pdf'), path = 'visualisations/N_delta/', height=12, width=15) 
  }
}


# Plot means vs observed values
for(c in n_country) {
  ggplot() +
    geom_point(data = P_df %>% filter(Country==c), aes(x=average_year, y=Observed, colour=Sector, pch=Sector)) +
    geom_errorbar(data = P_df %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax=upper_95, colour=Sector)) +
    geom_line(data = jags_Bspline_N_p_samps %>% filter(Country==c), aes(x=average_year, y=Mean, colour=Sector, lty=Source)) +
    geom_ribbon(data = jags_Bspline_N_p_samps %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
    geom_line(data = stan_Bspline_N_p_samps %>% filter(Country==c) , aes(x=average_year, y=Mean, colour=Sector, lty=Source)) +
    geom_ribbon(data = stan_Bspline_N_p_samps %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
    theme(title = element_text(size=10), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
    facet_wrap(~interaction(Method, Region), ncol = 5)
  ggsave(filename = paste0('compare_p_plot_',c,'_N_Bspline_JAGS_STAN.pdf'), path = 'visualisations/N_bspline/', height=12, width=15) 
  
}

# Plot means vs observed values
for(c in n_country) {
  for(s in n_sector) {
    print(c)
    print(c)
    ggplot() +
      geom_point(data = P_df %>% filter(Country==c & Sector==s), aes(x=average_year, y=Observed, pch=Sector)) +
      geom_errorbar(data = P_df %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax=upper_95)) +
      geom_line(data = original_P_samps %>% filter(Country==c & Sector==s), aes(x=average_year, y=Mean, colour=Source, lty=Source)) +
      geom_ribbon(data = original_P_samps %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax = upper_95,  lty=Source, fill=Source), alpha=0.2) +
      geom_line(data = jags_Bspline_N_p_samps %>% filter(Country==c & Sector==s), aes(x=average_year, y=Mean, colour=Source, lty=Source)) +
      geom_ribbon(data = jags_Bspline_N_p_samps %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax = upper_95, lty=Source, fill=Source), alpha=0.2) +
      # geom_line(data = stan_Bspline_N_p_samps %>% filter(Country==c & Sector==s) , aes(x=average_year, y=Mean, colour=Source, lty=Source)) +
      # geom_ribbon(data = stan_Bspline_N_p_samps %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax = upper_95, lty=Source, fill=Source), alpha=0.2) +
      theme(title = element_text(size=10), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
      facet_wrap(~interaction(Method, Region), ncol = 5)
    ggsave(filename = paste0('compare_p_plot_',c,'_',s,'_Bspline_JAGS_original.pdf'), path = 'visualisations/N_bspline/', height=12, width=15) 
  }
}


# Plot means vs observed values
for(c in n_country) {
  for(s in n_sector) {
    print(c)
    print(c)
    ggplot() +
      geom_point(data = P_df %>% filter(Country==c & Sector==s), aes(x=average_year, y=Observed, pch=Sector)) +
      geom_errorbar(data = P_df %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax=upper_95)) +
      geom_line(data = original_P_samps %>% filter(Country==c & Sector==s), aes(x=average_year, y=Mean, colour=Source, lty=Source)) +
      geom_ribbon(data = original_P_samps %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax = upper_95,  lty=Source, fill=Source), alpha=0.2) +
      geom_line(data = jags_MVN_p_samps %>% filter(Country==c & Sector==s), aes(x=average_year, y=Mean, colour=Source, lty=Source)) +
      geom_ribbon(data = jags_MVN_p_samps %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax = upper_95, lty=Source, fill=Source), alpha=0.2) +
      theme(title = element_text(size=10), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
      facet_wrap(~interaction(Method, Region), ncol = 5)
  ggsave(filename = paste0('compare_p_plot_',c,'_',s,'_MVN_JAGS_Original.pdf'), path = 'visualisations/MVN_delta/', height=12, width=15) 
  }
}

# Plot means vs observed values
for(c in n_country) {
  for(s in n_sector) {
    print(c)
    print(c)
    ggplot() +
      geom_point(data = P_df %>% filter(Country==c & Sector==s), aes(x=average_year, y=Observed, pch=Sector)) +
      geom_errorbar(data = P_df %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax=upper_95)) +
      geom_line(data = original_P_samps %>% filter(Country==c & Sector==s), aes(x=average_year, y=Mean, colour=Source, lty=Source)) +
      geom_ribbon(data = original_P_samps %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax = upper_95,  lty=Source, fill=Source), alpha=0.2) +
      geom_line(data = jags_Bspline_MVN_p_samps %>% filter(Country==c & Sector==s), aes(x=average_year, y=Mean, colour=Source, lty=Source)) +
      geom_ribbon(data = jags_Bspline_MVN_p_samps %>% filter(Country==c & Sector==s), aes(x=average_year, ymin=lower_95, ymax = upper_95, lty=Source, fill=Source), alpha=0.2) +
      theme(title = element_text(size=10), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
      facet_wrap(~interaction(Method, Region), ncol = 5)
    ggsave(filename = paste0('compare_p_plot_',c,'_',s,'_MVN_Bspline_JAGS_Original.pdf'), path = 'visualisations/MVN_bspline/', height=12, width=15) 
  }
}

