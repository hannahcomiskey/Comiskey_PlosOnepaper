library(tidyverse)
library(tidybayes)
library(mcmcplots)
library(bayesplot)
library(coda)
library(shinystan)
source('code/stan_utility.R')
source('code/load_functions.R')
source('code/DBDA2E-utilities.R')

vispath='visualisations/simulated_data/kstar/all_N/JAGS/beta_likelihood/'

# fit_sim model -----------------------------------------
fit_sim <- readRDS('results/JAGS_model_KenyaSim_allN_NCP_kstar_sum0_50years_beta.RDS')

fit_sim$model

plot(fit_sim)

model_samps <- fit_sim$BUGSoutput$sims.list

View(fit_sim$BUGSoutput$summary)

# Review spline coefficients
beta_k.mean <- apply(model_samps$beta.k, c(2,3,4), mean)
sum(beta_k.mean[1,1,1:13])

jags.fit.mcmc <- as.mcmc(fit_sim)

posterior <- as.data.frame(fit_sim$BUGSoutput$sims.list)

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(Year = all_years, index_year = 1:length(all_years))


# Plot traceplots -------------------------------------------------------------- 
mcmc_trace(jags.fit.mcmc, pars= c('alpha_cms[1,1]', 'alpha_cms[2,1]', 'alpha_cms[3,1]', 'alpha_cms[4,1]', 'alpha_cms[5,1]'))
mcmc_trace(jags.fit.mcmc, pars= c('alpha_pms[1,1]', 'alpha_pms[1,2]', 'alpha_pms[1,3]', 'alpha_pms[1,4]', 'alpha_pms[1,5]', 'alpha_pms[1,6]'))
mcmc_trace(jags.fit.mcmc, pars= c('alpha_pms[2,1]', 'alpha_pms[3,2]', 'alpha_pms[4,3]', 'alpha_pms[5,4]', 'alpha_pms[2,5]', 'alpha_pms[3,6]'))

mcmc_trace(jags.fit.mcmc, pars= c('sigma_delta')) #[1]', 'sigma_delta[2]', 'sigma_delta[3]')) #, 'sigma_delta[3]', 'sigma_delta[4]', 'sigma_delta[5]'))
mcmc_trace(jags.fit.mcmc, pars= c('sigma_alpha_pms[1]', 'sigma_alpha_pms[2]', 'sigma_alpha_pms[3]', 'sigma_alpha_pms[4]', 'sigma_alpha_pms[5]'))
mcmc_trace(jags.fit.mcmc, pars= c('sigma_alpha_cms[1]', 'sigma_alpha_cms[2]', 'sigma_alpha_cms[3]', 'sigma_alpha_cms[4]', 'sigma_alpha_cms[5]'))

mcmc_trace(jags.fit.mcmc, pars= c('delta.k[1,1,1]', 'delta.k[1,1,2]', 'delta.k[1,1,3]', 'delta.k[1,1,4]', 'delta.k[1,1,5]'))
mcmc_trace(jags.fit.mcmc, pars= c('delta.k[1,1,1]', 'delta.k[2,2,2]', 'delta.k[3,3,3]', 'delta.k[4,4,4]', 'delta.k[5,5,5]', 'delta.k[5,6,5]'))
mcmc_trace(jags.fit.mcmc, pars= c('delta.k[2,1,1]', 'delta.k[2,2,2]', 'delta.k[2,3,3]', 'delta.k[2,4,4]', 'delta.k[2,5,5]', 'delta.k[2,6,5]'))


# Check priors vs posteriors samples -------------------------------------------
beta_c_new <- t(as.matrix(beta_c_new))
colnames(beta_c_new) <- 1:5
beta_c_df <- as_tibble(beta_c_new) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  mutate(index_country = rep(1, each=1),
         Parameter = paste0("alpha_cms[",rep(1:5, 1),",",rep(1, each=5),"]"))

jags.fit.mcmc %>%
  spread_draws(`alpha_cms[5,1]`, `alpha_cms[4,1]`, `alpha_cms[3,1]`, `alpha_cms[2,1]`, `alpha_cms[1,1]`) %>%
  pivot_longer(cols = c(`alpha_cms[5,1]`, `alpha_cms[4,1]`, `alpha_cms[3,1]`, `alpha_cms[2,1]`, `alpha_cms[1,1]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  geom_vline(data=beta_c_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = "betac_plot_Kenyasim_50years.pdf", path = vispath, height=12, width=15) 


alpha_sim <- t(alpha_sim)
colnames(alpha_sim) <- 1:5
alpha_df <- as_tibble(alpha_sim) %>%
  mutate(index_subnat = 1:P) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  rowwise() %>%
  mutate(Parameter = paste0("alpha_pms[", index_method,",",index_subnat,"]"))


jags.fit.mcmc %>%
  spread_draws(`alpha_pms[1,1]`, `alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[1,6]`,
               `alpha_pms[1,2]`, `alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[2,6]`,
               `alpha_pms[1,3]`, `alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[3,6]`,
               `alpha_pms[1,4]`, `alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[4,6]`,
               `alpha_pms[1,5]`, `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[5,6]`
  ) %>%
  pivot_longer(cols = c(`alpha_pms[1,1]`, `alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[1,6]`,
                        `alpha_pms[1,2]`, `alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[2,6]`,
                        `alpha_pms[1,3]`, `alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[3,6]`,
                        `alpha_pms[1,4]`, `alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[4,6]`,
                        `alpha_pms[1,5]`, `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[5,6]`
                        ), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  geom_vline(data=alpha_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = "alpha_pms_plot_Kenyasim_50years.pdf", path = vispath, height=12, width=15) 

# Posterior predictive check
y_sim <- fit_sim$BUGSoutput$sims.list$y.sim
y_sim_mean <- apply(y_sim, 2, mean)
y_obs <- as.vector(unlist(logit.data[,c("logit.Public")]))

pdf(paste0(vispath, "posterior_predictive_observed_Kenyasim_50years.pdf")) 
plot(y_sim_mean, y_obs)
abline(a=0, b=1, col='red')
dev.off()

pdf(paste0(vispath,"posterior_predictive_check_Kenyasim_50years.pdf"))
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

P_samps <- model_samps$P
dim(P_samps)
P_samps.mean <- apply(P_samps, c(2,3,4,5), mean)
P_samps.mean <- plyr::adply(P_samps.mean, .margins=c(2,3,4))
colnames(P_samps.mean) <- c('index_method', 'index_subnat', 'index_year', 'Public', 'Private') 
P_samps.mean <- P_samps.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')

P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
dim(P_samps.quantile)
colnames(P_samps.quantile) <- c('index_sector', 'index_method', 'index_subnat', 'index_year', 'lower_95', 'upper_95') 
P_samps.quantile <- P_samps.quantile %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  left_join(sector_index_table) 

P_samps_df <- left_join(P_samps.mean, P_samps.quantile)
# saveRDS(P_samps_df, file='results/JAGS_model_Kenya_MVN_NCP_kstar_dsigma_10years_Psamps_df.RDS')

# Get observed data 
P_df<- P_sim_df_sample %>%
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

# Plot means vs observed values
ggplot() +
  geom_point(data = P_df, aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_samps_df  , aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df , aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  facet_wrap(~interaction(Method, index_subnat), ncol=5)
ggsave(filename = "p_plot_Kenyasim_sum0_50years.pdf", path = vispath, height=12, width=15) 

alpha_samps <- model_samps$alpha_pms
dim(alpha_samps)
alpha_mean <- apply(alpha_samps, c(2,3), mean)
colnames(alpha_mean) <- 1:P

alpha_mean <- as_tibble(alpha_mean) %>%
  mutate(index_method = 1:M) %>%
  pivot_longer(cols=`1`:`6`, names_to = 'index_subnat', values_to = 'alpha') %>%
  mutate(invlogit.alpha = exp(alpha)/(1+exp(alpha))) %>%
  left_join(method_index_table)

ggplot() +
  geom_point(data = P_df , aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_samps_df  , aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df , aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  geom_hline(data=alpha_mean, aes(yintercept = invlogit.alpha)) +
  facet_wrap(~interaction(Method, index_subnat), ncol = 5)
ggsave(filename = "p_plot_Kenyasim_10years_withAlpha_dsigma.pdf", path = vispath, height=12, width=15) 
