library(tidyverse)
library(tidybayes)
library(mcmcplots)
library(bayesplot)
library(coda)
library(shinystan)
library(R2jags)
library(rjags)
source('code/load_functions.R')
source('code/DBDA2E-utilities.R')

vispath='visualisations/simulated_data/MVN_bspline/JAGS/'

# fit_sim model -----------------------------------------
fit_sim <- readRDS('results/JAGS/JAGS_model_KenyaSim_MVN_Bspline_10years.RDS')

fit_sim$model

plot(fit_sim)

model_samps <- fit_sim$BUGSoutput$sims.list

View(fit_sim$BUGSoutput$summary)

# Review spline coefficients
a.mean <- apply(model_samps$a, c(2,3,4), mean)
sum(a.mean[1,1,1:13])

jags.fit.mcmc <- as.mcmc(fit_sim)

posterior <- as.data.frame(fit_sim$BUGSoutput$sims.list)

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(Year = all_years, index_year = 1:length(all_years))

# Plot traceplots -------------------------------------------------------------- 
mcmc_trace(jags.fit.mcmc, pars= c('beta_c[1,1]', 'beta_c[2,1]', 'beta_c[3,1]', 'beta_c[4,1]', 'beta_c[5,1]'))

mcmc_trace(jags.fit.mcmc, pars= c('err_a0[1,1]', 'err_a0[2,1]', 'err_a0[3,1]', 'err_a0[4,1]', 'err_a0[5,1]'))
mcmc_trace(jags.fit.mcmc, pars= c('err_a0[1,1]', 'err_a0[1,2]', 'err_a0[1,3]', 'err_a0[1,4]', 'err_a0[1,5]'))

mcmc_trace(jags.fit.mcmc, pars= c('err_a[1,1,1]', 'err_a[2,1,2]', 'err_a[3,1,3]', 'err_a[4,1,4]', 'err_a[5,1,5]'))
mcmc_trace(jags.fit.mcmc, pars= c('err_a[1,1,1]', 'err_a[2,2,2]', 'err_a[3,3,3]', 'err_a[4,4,4]', 'err_a[5,5,5]'))

mcmc_trace(jags.fit.mcmc, pars= c('sigma_y[1]', 'sigma_y[2]', 'sigma_y[3]', 'sigma_y[4]', 'sigma_y[5]'))

mcmc_trace(jags.fit.mcmc, pars= c('a[1,1,1]', 'a[1,1,2]', 'a[1,1,3]', 'a[1,1,4]', 'a[1,1,5]'))

# Check priors vs posteriors samples -------------------------------------------
beta_c_new <- t(as.matrix(beta_c_new))
colnames(beta_c_new) <- 1:5
beta_c_df <- as_tibble(beta_c_new) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  mutate(index_country = rep(1, each=1),
         Parameter = paste0("beta_c[",1:5,",1]"))

jags.fit.mcmc %>%
  spread_draws(`beta_c[1,1]`, `beta_c[2,1]`, `beta_c[3,1]`, `beta_c[4,1]`, `beta_c[5,1]`) %>%
  pivot_longer(cols = c(`beta_c[1,1]`, `beta_c[2,1]`, `beta_c[3,1]`, `beta_c[4,1]`, `beta_c[5,1]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  geom_vline(data=beta_c_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = "betac_plot_Kenyasim_10years_mvn_bspline.pdf", path = vispath, height=12, width=15) 

alpha_sim <- t(alpha_sim)
colnames(alpha_sim) <- 1:5
alpha_df <- as_tibble(alpha_sim) %>%
  mutate(index_subnat = 1:P) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  rowwise() %>%
  mutate(Parameter = paste0("a0[", index_method,",",index_subnat,"]"))


jags.fit.mcmc %>%
  spread_draws(`a0[1,1]`, `a0[2,1]`, `a0[3,1]`, `a0[4,1]`, `a0[5,1]`, `a0[1,6]`,
               `a0[1,2]`, `a0[2,2]`, `a0[3,2]`, `a0[4,2]`, `a0[5,2]`, `a0[2,6]`,
               `a0[1,3]`, `a0[2,3]`, `a0[3,3]`, `a0[4,3]`, `a0[5,3]`, `a0[3,6]`,
               `a0[1,4]`, `a0[2,4]`, `a0[3,4]`, `a0[4,4]`, `a0[5,4]`, `a0[4,6]`,
               `a0[1,5]`, `a0[2,5]`, `a0[3,5]`, `a0[4,5]`, `a0[5,5]`, `a0[5,6]`
  ) %>%
  pivot_longer(cols = c(`a0[1,1]`, `a0[2,1]`, `a0[3,1]`, `a0[4,1]`, `a0[5,1]`, `a0[1,6]`,
                        `a0[1,2]`, `a0[2,2]`, `a0[3,2]`, `a0[4,2]`, `a0[5,2]`, `a0[2,6]`,
                        `a0[1,3]`, `a0[2,3]`, `a0[3,3]`, `a0[4,3]`, `a0[5,3]`, `a0[3,6]`,
                        `a0[1,4]`, `a0[2,4]`, `a0[3,4]`, `a0[4,4]`, `a0[5,4]`, `a0[4,6]`,
                        `a0[1,5]`, `a0[2,5]`, `a0[3,5]`, `a0[4,5]`, `a0[5,5]`, `a0[5,6]`
                        ), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  geom_vline(data=alpha_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = "a0_plot_Kenyasim_10years_mvn_bspline.pdf", path = vispath, height=12, width=15) 

# Posterior predictive check
y_sim <- fit_sim$BUGSoutput$sims.list$Y_sim
y_sim_mean <- apply(y_sim, 2, mean)
y_obs <- as.vector(unlist(logit.data[,c("logit.Public")]))

pdf(paste0(vispath, "posterior_predictive_observed_Kenyasim_10years_mvn_bspline.pdf")) 
plot(y_sim_mean, y_obs)
abline(a=0, b=1, col='red')
dev.off()

pdf(paste0(vispath,"posterior_predictive_check_Kenyasim_10years_mvn_bspline.pdf"))
# plot the observed data hist
hist(y_obs, freq = FALSE, breaks = 50) # observed data 
n_samples = nrow(y_sim) 
for (r in 1:n_samples){ 
  lines(density(y_sim[r,]), col = rgb(135, 206, 235, max= 255, alpha=25)) 
}
dev.off()

# Get P estimates --------------------------------------------------------------

P_samps <- model_samps$P
dim(P_samps)
P_samps.mean <- apply(P_samps, c(2,3,4,5), mean)
P_samps.mean <- plyr::adply(P_samps.mean, .margins=c(1,2,4))
colnames(P_samps.mean) <- c('index_subnat', 'index_method', 'index_year', 'Public', 'Private') 
P_samps.mean <- P_samps.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')

P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
dim(P_samps.quantile)
colnames(P_samps.quantile) <- c('index_subnat', 'index_method', 'index_sector', 'index_year', 'lower_95', 'upper_95') 
P_samps.quantile <- P_samps.quantile %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  left_join(sector_index_table) 

P_samps_df <- left_join(P_samps.mean, P_samps.quantile)
# saveRDS(P_samps_df, file='results/JAGS_model_Kenya_Bspline_10years_Psamps_df.RDS')

# Get observed data 
P_df<- P_sim_df_sample %>%
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

# Plot means vs observed values
ggplot() +
  geom_point(data = P_df, aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_samps_df , aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df, aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  facet_wrap(~interaction(Method, index_subnat), ncol=5)
ggsave(filename = "p_plot_Kenyasim_sum0_10years_mvn_bspline.pdf", path = vispath, height=12, width=15) 

alpha_samps <- model_samps$a0
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
  geom_line(data = P_samps_df, aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df, aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  geom_hline(data=alpha_mean, aes(yintercept = invlogit.alpha)) +
  facet_wrap(~interaction(Method, index_subnat), ncol = 5)
ggsave(filename = "p_plot_Kenyasim_10years_withAlpha_bspline.pdf", path = vispath, height=12, width=15) 
