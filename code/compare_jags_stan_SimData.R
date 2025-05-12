library(rjags)
library(R2jags)
library(rstan)
library(tidyverse)
library(tidybayes)
library(parallel)
source('code/load_functions.R')

# Source simulated data --------------------------------------
load("data/simulated_data/simulated_data_all_N_kstar_Kenya_10years.RData")

# Get logit of parameters and variance -----------------------------------------
mydata <- P_sim_df_sample[,c("Public")] %>%
  mutate(Public.SE = 0.1) %>%
  rowwise() %>%
  mutate(Public = ifelse(Public < 0.0001 , 0, Public))

logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)))

# # testing splines ------------------------------------------------------------
all_years <-  -5:55 # -5:20 #
B <- bs_bbase_precise(all_years)
Bik <- B$B.ik
K <-dim(Bik)[2]
H = K-1
kstar = B$Kstar

year_index_table <- tibble(Year = all_years, index_year = 1:length(all_years))

P_sim_df_sample <- P_sim_df_sample %>% 
  rename(Year = index_year) %>%
  mutate_if(is.character, as.numeric) %>%
  left_join(year_index_table)

# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- matchcountry
n_all_years <- length(all_years)
M_count = 5

# Read in results --------------------------------------------------------------
STAN_sim <- readRDS('results/STAN_model_allN_NCP_kstar_sum0_KenyaSim_10years_dsigma.RDS')
stan_samps <- rstan::extract(STAN_sim)

JAGS_sim <- readRDS('results/JAGS_model_KenyaSim_allN_NCP_kstar_sum0_10years_dsigma.RDS')
jags_samps <- JAGS_sim$BUGSoutput$sims.list

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(Year = all_years, index_year = 1:length(all_years))

# Compare traceplots of sigma_delta --------------------------------------------
chain <- rep(1:3, each=2000)
sigma_stan <- stan_samps$sigma_delta
sigma_jags <- jags_samps$sigma_delta

dataframe <- tibble(STAN=sigma_stan, JAGS=as.vector(sigma_jags), Chain = as.factor(chain), iteration = rep(1:2000, 3))

dataframe_long <- dataframe %>% pivot_longer(cols=c(JAGS, STAN), names_to = 'Source', values_to = 'Sample')

ggplot(data = dataframe_long) +
  geom_line(aes(x=iteration, y=Sample, colour = Chain), alpha=0.7) +
  theme(legend.position = 'bottom') +
  facet_wrap(~Source)
ggsave(filename = "sigma_delta_JAGS_STAN_Kenyasim_10years.pdf", path = 'visualisations/simulated_data/kstar/all_N/', height=12, width=15) 

# Compare traceplots of sigma_alpha --------------------------------------------
sigma_stan <- stan_samps$sigma_alpha
colnames(sigma_stan) <- n_method
sigma_stan <- as_tibble(sigma_stan) %>% 
  mutate(Source = 'STAN', Chain = as.factor(chain), iteration = rep(1:2000, 3)) %>%
  pivot_longer(cols=c(`Female Sterilization`, `Implants`, `Injectables`,   `IUD`, `OC Pills`), names_to = 'Method', values_to = 'Sample') 

sigma_jags <- jags_samps$sigma_alpha_pms
colnames(sigma_jags) <- n_method
sigma_jags <- as_tibble(sigma_jags) %>% 
  mutate(Source = 'JAGS', Chain = as.factor(chain), iteration = rep(1:2000, 3)) %>%
  pivot_longer(cols=c(`Female Sterilization`, `Implants`, `Injectables`,   `IUD`, `OC Pills`), names_to = 'Method', values_to = 'Sample')

dataframe <- tibble(rbind(sigma_jags, sigma_stan))

ggplot(data = dataframe) +
  geom_line(aes(x=iteration, y=Sample, colour = Chain), alpha=0.7) +
  theme(legend.position = 'bottom') +
  facet_wrap(~ interaction(Source, Method), ncol=2, scales='free_y')
ggsave(filename = "sigma_alpha_pms_JAGS_STAN_Kenyasim_10years.pdf", path = 'visualisations/simulated_data/kstar/all_N/', height=12, width=15) 

# Compare traceplots of alpha_cms --------------------------------------------
alpha_cms_stan <- stan_samps$beta_c
alpha_cms_stan <- plyr::adply(alpha_cms_stan, c(2))
alpha_cms_stan <- alpha_cms_stan[,2:6]
colnames(alpha_cms_stan) <- n_method
alpha_cms_stan <- as_tibble(alpha_cms_stan) %>% 
  mutate(Source = 'STAN', Chain = as.factor(chain), iteration = rep(1:2000, 3)) %>%
  pivot_longer(cols=c(`Female Sterilization`, `Implants`, `Injectables`,   `IUD`, `OC Pills`), names_to = 'Method', values_to = 'Sample') 

alpha_cms_jags <- jags_samps$alpha_cms
alpha_cms_jags <- plyr::adply(alpha_cms_jags, c(3))
alpha_cms_jags <- alpha_cms_jags[,2:6]
colnames(alpha_cms_jags) <- n_method
alpha_cms_jags <- as_tibble(alpha_cms_jags) %>% 
  mutate(Source = 'JAGS', Chain = as.factor(chain), iteration = rep(1:2000, 3)) %>%
  pivot_longer(cols=c(`Female Sterilization`, `Implants`, `Injectables`,   `IUD`, `OC Pills`), names_to = 'Method', values_to = 'Sample')

dataframe <- tibble(rbind(alpha_cms_jags, alpha_cms_stan))

ggplot(data = dataframe) +
  geom_line(aes(x=iteration, y=Sample, colour = Chain), alpha=0.7) +
  theme(legend.position = 'bottom') +
  facet_wrap(~ interaction(Source, Method), ncol=2, scales='free_y')
ggsave(filename = "alpha_cms_JAGS_STAN_Kenyasim_10years.pdf", path = 'visualisations/simulated_data/kstar/all_N/', height=12, width=15) 


# P ----------------------------------------------------------------------------
P_jags <- jags_samps$P
dim(P_jags)
P_jags.mean <- apply(P_jags, c(2,3,4,5), mean)
P_jags.mean <- plyr::adply(P_jags.mean, .margins=c(2,3,4))
colnames(P_jags.mean) <- c('index_method', 'index_subnat', 'index_year', 'Public', 'Private') 
P_jags.mean <- P_jags.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')

P_jags.quantile <- apply(P_jags, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
P_jags.quantile <- plyr::adply(P_jags.quantile, .margins=c(2,3,4,5))
dim(P_jags.quantile)
colnames(P_jags.quantile) <- c('index_sector', 'index_method', 'index_subnat', 'index_year', 'lower_95', 'upper_95') 
P_jags.quantile <- P_jags.quantile %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  left_join(sector_index_table) 

P_jags_df <- left_join(P_jags.mean, P_jags.quantile) %>%
  mutate(Source = 'JAGS')

P_stan <- stan_samps$P
dim(P_stan)
P_stan.mean <- apply(P_stan, c(2,3,4,5), mean)
P_stan.mean <- plyr::adply(P_stan.mean, .margins=c(1,2,4))
colnames(P_stan.mean) <- c('index_method', 'index_subnat', 'index_year', 'Public', 'Private') 
P_stan.mean <- P_stan.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')

P_stan <- stan_samps$P
dim(P_stan)
P_stan.quantile <- apply(P_stan, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
P_stan.quantile <- plyr::adply(P_stan.quantile, .margins=c(2,3,4,5))
colnames(P_stan.quantile) <- c('index_method', 'index_subnat','index_sector', 'index_year', 'lower_95', 'upper_95') 
P_stan.quantile <- P_stan.quantile %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  left_join(sector_index_table) 

P_stan_df <- left_join(P_stan.mean, P_stan.quantile) %>%
  mutate(Source = 'STAN')

dataframe <- tibble(rbind(P_stan_df, P_jags_df))

# Get observed data 
P_df<- P_sim_df_sample %>%
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

# Plot means vs observed values
ggplot() +
  geom_point(data = P_df, aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = dataframe %>% filter(index_year <20), aes(x=index_year, y=Mean, colour=Sector, lty=Source)) +
  # geom_ribbon(data = dataframe %>% filter(index_year <20) , aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  facet_wrap(~interaction(Method, index_subnat), ncol=5)
ggsave(filename = "p_plot_JAGS_STAN_Kenyasim_10years.pdf", path = 'visualisations/simulated_data/kstar/all_N/', height=17, width=20) 

