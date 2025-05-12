library(tidyverse)
library(tidybayes)
library(mcmcplots)
library(bayesplot)
library(BayesTools)
library(coda)
library(shinystan)

# Source code --------------------------------------
source('code/stan_utility.R')
source('code/load_functions.R')
source('code/DBDA2E-utilities.R')
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

# fit_sim model -----------------------------------------
vispath = 'visualisations/N_delta/'

fit_sim <- readRDS('results/JAGS_mod_N_delta_all_parallel.RDS')

fit_sim$model

plot(fit_sim)

print(fit_sim)

View(fit_sim$BUGSoutput$summary)

model_samps <- fit_sim$BUGSoutput$sims.list

# Review spline coefficients
beta_k.mean <- apply(model_samps$beta.k, c(2,3,4), mean)
sum(beta_k.mean[1,1,1:13])

jags.fit.mcmc <- as.mcmc(fit_sim)

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))
subnat_index_table <- FP_source_data_wide %>% select(Country, Region, index_country, index_subnat) %>% distinct()

# Plot ACF ---------------------------------------------------------------------
mcmc_acf(jags.fit.mcmc, pars='sigma_delta')
mcmc_acf(jags.fit.mcmc, pars=c('sigma_alpha_cms[1]'))
mcmc_acf(jags.fit.mcmc, pars=c('sigma_alpha_pms[1]'))

# Plot traceplots -------------------------------------------------------------- 
mcmc_trace(jags.fit.mcmc, pars= c('alpha_cms[1,1]', 'alpha_cms[2,1]', 'alpha_cms[3,1]', 'alpha_cms[4,1]', 'alpha_cms[5,1]',
                                  'alpha_cms[1,2]', 'alpha_cms[2,2]', 'alpha_cms[3,2]', 'alpha_cms[4,2]', 'alpha_cms[5,2]',
                                  'alpha_cms[1,3]', 'alpha_cms[2,3]', 'alpha_cms[3,3]', 'alpha_cms[4,3]', 'alpha_cms[5,3]',
                                  'alpha_cms[1,4]', 'alpha_cms[2,4]', 'alpha_cms[3,4]', 'alpha_cms[4,4]', 'alpha_cms[5,4]'))
                                  # 'alpha_cms[1,5]', 'alpha_cms[2,5]', 'alpha_cms[3,5]', 'alpha_cms[4,5]', 'alpha_cms[5,5]',
                                  # 'alpha_cms[1,6]', 'alpha_cms[2,6]', 'alpha_cms[3,6]', 'alpha_cms[4,6]', 'alpha_cms[5,6]',
                                  # 'alpha_cms[1,7]', 'alpha_cms[2,7]', 'alpha_cms[3,7]', 'alpha_cms[4,7]', 'alpha_cms[5,7]'))

mcmc_trace(jags.fit.mcmc, pars= c('alpha_pms[1,1]', 'alpha_pms[1,2]', 'alpha_pms[1,3]', 'alpha_pms[1,4]', 'alpha_pms[1,5]', 'alpha_pms[1,6]'))
mcmc_trace(jags.fit.mcmc, pars= c('alpha_pms[2,1]', 'alpha_pms[3,2]', 'alpha_pms[4,3]', 'alpha_pms[5,4]', 'alpha_pms[2,5]', 'alpha_pms[3,6]'))

mcmc_trace(jags.fit.mcmc, pars= c('sigma_delta')) #[1]', 'sigma_delta[2]', 'sigma_delta[3]', 'sigma_delta[4]', 'sigma_delta[5]'))

mcmc_trace(jags.fit.mcmc, pars= c('delta.k[1,1,1]', 'delta.k[1,1,2]', 'delta.k[1,1,3]', 'delta.k[1,1,4]', 'delta.k[1,1,5]'))
mcmc_trace(jags.fit.mcmc, pars= c('delta.k[1,1,1]', 'delta.k[2,2,2]', 'delta.k[3,3,3]', 'delta.k[4,4,4]', 'delta.k[5,5,5]', 'delta.k[5,6,5]'))

mcmc_trace(jags.fit.mcmc, pars= c('delta.k[2,1,1]', 'delta.k[2,2,2]', 'delta.k[2,3,3]', 'delta.k[2,4,4]', 'delta.k[2,5,5]', 'delta.k[2,6,5]'))

mcmc_trace(jags.fit.mcmc, pars= c('sigma_alpha_pms[1]', 'sigma_alpha_pms[2]', 'sigma_alpha_pms[3]', 'sigma_alpha_pms[4]', 'sigma_alpha_pms[5]'))

# Check priors vs posteriors samples -------------------------------------------

jags.fit.mcmc %>%
  spread_draws(`alpha_cms[5,1]`, `alpha_cms[4,1]`, `alpha_cms[3,1]`, `alpha_cms[2,1]`, `alpha_cms[1,1]`,
               `alpha_cms[5,2]`, `alpha_cms[4,2]`, `alpha_cms[3,2]`, `alpha_cms[2,2]`, `alpha_cms[1,2]`,
               `alpha_cms[5,3]`, `alpha_cms[4,3]`, `alpha_cms[3,3]`, `alpha_cms[2,3]`, `alpha_cms[1,3]`,
               `alpha_cms[5,4]`, `alpha_cms[4,4]`, `alpha_cms[3,4]`, `alpha_cms[2,4]`, `alpha_cms[1,4]`) %>%
               # `alpha_cms[5,5]`, `alpha_cms[4,5]`, `alpha_cms[3,5]`, `alpha_cms[2,5]`, `alpha_cms[1,5]`,
               # `alpha_cms[5,6]`, `alpha_cms[4,6]`, `alpha_cms[3,6]`, `alpha_cms[2,6]`, `alpha_cms[1,6]`,
               # `alpha_cms[5,7]`, `alpha_cms[4,7]`, `alpha_cms[3,7]`, `alpha_cms[2,7]`, `alpha_cms[1,7]`) %>%
               # `alpha_cms[5,8]`, `alpha_cms[4,8]`, `alpha_cms[3,8]`, `alpha_cms[4,8]`, `alpha_cms[1,8]`) %>%
  pivot_longer(cols = c(`alpha_cms[5,1]`, `alpha_cms[4,1]`, `alpha_cms[3,1]`, `alpha_cms[2,1]`, `alpha_cms[1,1]`,
                        `alpha_cms[5,2]`, `alpha_cms[4,2]`, `alpha_cms[3,2]`, `alpha_cms[2,2]`, `alpha_cms[1,2]`,
                        `alpha_cms[5,3]`, `alpha_cms[4,3]`, `alpha_cms[3,3]`, `alpha_cms[2,3]`, `alpha_cms[1,3]`,
                        `alpha_cms[5,4]`, `alpha_cms[4,4]`, `alpha_cms[3,4]`, `alpha_cms[2,4]`, `alpha_cms[1,4]`), names_to='Parameter', values_to = 'sample') %>%
                        # `alpha_cms[5,5]`, `alpha_cms[4,5]`, `alpha_cms[3,5]`, `alpha_cms[2,5]`, `alpha_cms[1,5]`
                        # `alpha_cms[5,6]`, `alpha_cms[4,6]`, `alpha_cms[3,6]`, `alpha_cms[2,6]`, `alpha_cms[1,6]`,
                        # `alpha_cms[5,7]`, `alpha_cms[4,7]`, `alpha_cms[3,7]`, `alpha_cms[2,7]`, `alpha_cms[1,7]`
                        # `alpha_cms[5,8]`, `alpha_cms[4,8]`, `alpha_cms[3,8]`, `alpha_cms[4,8]`, `alpha_cms[1,8]`
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  theme(title = element_text(size=20), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 20), legend.text = element_text(size = 20))+
  facet_wrap(~Parameter, ncol=5)
ggsave(filename = "/parameter_plots/betac_plot_allN_NCP.pdf", path = vispath, height=12, width=15) 


jags.fit.mcmc %>%
  spread_draws(`alpha_pms[1,1]`, `alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[1,6]`, `alpha_pms[1,7]`, `alpha_pms[1,8]`, `alpha_pms[1,9]`, `alpha_pms[1,10]`, 
               `alpha_pms[1,2]`, `alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[2,6]`, `alpha_pms[2,7]`, `alpha_pms[2,8]`, `alpha_pms[2,9]`, `alpha_pms[1,10]`,
               `alpha_pms[1,3]`, `alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[3,6]`, `alpha_pms[3,7]`, `alpha_pms[3,8]`, `alpha_pms[3,9]`, `alpha_pms[1,10]`,
               `alpha_pms[1,4]`, `alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[4,6]`, `alpha_pms[4,7]`, `alpha_pms[4,8]`, `alpha_pms[4,9]`, `alpha_pms[1,10]`,
               `alpha_pms[1,5]`, `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[5,6]`, `alpha_pms[5,7]`, `alpha_pms[5,8]`, `alpha_pms[5,9]`, `alpha_pms[1,10]`) %>%
  pivot_longer(cols = c(`alpha_pms[1,1]`, `alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[1,6]`, `alpha_pms[1,7]`, `alpha_pms[1,8]`, `alpha_pms[1,9]`, `alpha_pms[1,10]`, 
                        `alpha_pms[1,2]`, `alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[2,6]`, `alpha_pms[2,7]`, `alpha_pms[2,8]`, `alpha_pms[2,9]`, `alpha_pms[1,10]`,
                        `alpha_pms[1,3]`, `alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[3,6]`, `alpha_pms[3,7]`, `alpha_pms[3,8]`, `alpha_pms[3,9]`, `alpha_pms[1,10]`,
                        `alpha_pms[1,4]`, `alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[4,6]`, `alpha_pms[4,7]`, `alpha_pms[4,8]`, `alpha_pms[4,9]`, `alpha_pms[1,10]`,
                        `alpha_pms[1,5]`, `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[5,6]`, `alpha_pms[5,7]`, `alpha_pms[5,8]`, `alpha_pms[5,9]`, `alpha_pms[1,10]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  theme(title = element_text(size=20), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 20), legend.text = element_text(size = 20))+
  facet_wrap(~Parameter, ncol=5)
ggsave(filename = "/parameter_plots/alpha_pms_plot_allN_NCP.pdf", path = vispath, height=12, width=15) 

# Posterior predictive check
y_sim <- fit_sim$BUGSoutput$sims.list$y.sim
y_sim_mean <- apply(y_sim, 2, mean)
y_obs <- as.vector(unlist(logit.data[,c("logit.Public")]))

pdf(paste0(vispath, "/parameter_plots/posterior_predictive_observed_4country_allN_NCP.pdf"))
plot(y_sim_mean, y_obs)
abline(a=0, b=1, col='red')
dev.off()

pdf(paste0(vispath, "/parameter_plots/posterior_predictive_check_4country_allN_NCP.pdf"))
# plot the observed data hist
hist(y_obs, freq = FALSE, breaks = 50) # observed data 
# number of samples to simulate 
n_samples = nrow(y_sim) 
#simulate samples 
for (r in 1:n_samples){ 
  # add plot of simulated sample to histogram 
  lines(density(y_sim[r,]), col = rgb(135, 206, 235, max= 255, alpha=25)) 
}
dev.off()

# Get P estimates --------------------------------------------------------------
# P_samps <- fit_sim$BUGSoutput$sims.list$P
# dim(P_samps)
# P_samps.mean <- apply(P_samps, c(2,3,4,5), mean)
# P_samps.mean <- plyr::adply(P_samps.mean, .margins=c(2,3,4))
# colnames(P_samps.mean) <- c('index_method', 'index_subnat', 'index_year', 'Public', 'Private') 
# P_samps.mean <- P_samps.mean %>% 
#   mutate(across(everything(), as.numeric)) %>%
#   left_join(method_index_table) %>%
#   pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')
# 
# P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
# P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
# dim(P_samps.quantile)
# colnames(P_samps.quantile) <- c('index_sector', 'index_method', 'index_subnat', 'index_year', 'lower_95', 'upper_95') 
# P_samps.quantile <- P_samps.quantile %>% 
#   mutate(across(everything(), as.numeric)) %>%
#   left_join(method_index_table) %>%
#   left_join(sector_index_table) %>%
#   left_join(index_country_subnat_tbl) %>%
#   left_join(year_index_table)
# 
# P_samps_df <- left_join(P_samps.mean, P_samps.quantile)
# saveRDS(P_samps_df, file='results/JAGS/P_samps/JAGS_allN_Psamps_df.RDS')

P_samps_df <- readRDS('results/JAGS/P_samps/N_delta/JAGS_N_delta_Psamps_df.RDS')

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

# Plot means vs observed values
for(c in n_country) {
  ggplot() +
    geom_point(data = P_df %>% filter(Country==c), aes(x=average_year, y=Observed, colour=Sector, pch=Sector)) +
    geom_errorbar(data = P_df %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax=upper_95, colour=Sector)) +
    geom_line(data = P_samps_df %>% filter(Country==c), aes(x=average_year, y=Mean, colour=Sector, lty=Sector)) +
    geom_ribbon(data = P_samps_df %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
    theme(title = element_text(size=20), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=8), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
    theme(legend.position = "bottom", legend.title = element_text(size = 20), legend.text = element_text(size = 20))+
    facet_wrap(~interaction(Region, Method), ncol=5)
  ggsave(filename = paste0("p_plot_",c,".pdf"), path = paste0(vispath,'/country_plots/'), height=12, width=15) 
  
}

alpha_samps <- model_samps$alpha_pms
dim(alpha_samps)
alpha_mean <- apply(alpha_samps, c(2,3), mean)
colnames(alpha_mean) <- n_subnat

alpha_mean <- as_tibble(alpha_mean) %>%
  mutate(index_method = 1:length(n_method)) %>%
  pivot_longer(cols=all_of(n_subnat), names_to = 'Region', values_to = 'alpha') %>%
  mutate(invlogit.alpha = exp(alpha)/(1+exp(alpha))) %>%
  left_join(method_index_table) %>%
  left_join(index_country_subnat_tbl)

# Plot means vs observed values
for(c in n_country) {
  ggplot() +
    geom_point(data = P_df %>% filter(Country==c), aes(x=average_year, y=Observed, colour=Sector, pch=Sector)) +
    geom_errorbar(data = P_df %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax=upper_95, colour=Sector)) +
    geom_line(data = P_samps_df %>% filter(Country==c), aes(x=average_year, y=Mean, colour=Sector, lty=Sector)) +
    geom_ribbon(data = P_samps_df %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
    geom_hline(data=alpha_mean %>% filter(Country==c), aes(yintercept = invlogit.alpha)) +
    theme(title = element_text(size=20), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
    theme(legend.position = "bottom", legend.title = element_text(size = 20), legend.text = element_text(size = 20))+
    facet_wrap(~interaction(Method, Region), ncol=5)
  ggsave(filename = paste0("p_plot_",c,"_withAlpha.pdf"), path =  paste0(vispath,'/country_plots/'), height=12, width=15) 
}

