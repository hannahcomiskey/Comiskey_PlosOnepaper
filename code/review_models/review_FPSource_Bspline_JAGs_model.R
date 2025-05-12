library(tidyverse)
library(tidybayes)
library(mcmcplots)
library(bayesplot)
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
vispath = 'visualisations/country_data/bspline/JAGS/TG_spline/'

# fit_sim model -----------------------------------------
fit_sim <- readRDS('results/JAGS/JAGS_mod_Bspline_N_TG_all.RDS')

fit_sim$model

plot(fit_sim)

# model_samps <- fit_sim$BUGSoutput$sims.list

jags.fit.mcmc <- as.mcmc(fit_sim)

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))

# Plot traceplots -------------------------------------------------------------- 
mcmc_trace(jags.fit.mcmc, pars= c('a_xi', 'c_xi', 'kappa_sq_B'))

mcmc_trace(jags.fit.mcmc, pars= c('xi_sq[1]', 'xi_sq[2]', 'xi_sq[3]', 'xi_sq[4]', 'xi_sq[5]'))

mcmc_trace(jags.fit.mcmc, pars= c('kappa_sq[1]', 'kappa_sq[2]', 'kappa_sq[3]', 'kappa_sq[4]', 'kappa_sq[5]'))

mcmc_trace(jags.fit.mcmc, pars= c('beta_c[1,1]', 'beta_c[2,1]', 'beta_c[3,1]', 'beta_c[4,1]', 'beta_c[5,1]',
                                  'beta_c[1,2]', 'beta_c[2,2]', 'beta_c[3,2]', 'beta_c[4,2]', 'beta_c[5,2]',
                                  'beta_c[1,3]', 'beta_c[2,3]', 'beta_c[3,3]', 'beta_c[4,3]', 'beta_c[5,3]',
                                  'beta_c[1,4]', 'beta_c[2,4]', 'beta_c[3,4]', 'beta_c[4,4]', 'beta_c[5,4]'))

mcmc_trace(jags.fit.mcmc, pars= c('sigma_a0[1]', 'sigma_a0[2]', 'sigma_a0[3]', 'sigma_a0[4]', 'sigma_a0[5]'))

mcmc_trace(jags.fit.mcmc, pars= c('sigma_a[1]', 'sigma_a[2]', 'sigma_a[3]', 'sigma_a[4]', 'sigma_a[5]'))

mcmc_trace(jags.fit.mcmc, pars= c('a[1,1,1]', 'a[1,1,2]', 'a[1,1,3]', 'a[1,1,4]', 'a[1,1,5]'))

# Check priors vs posteriors samples -------------------------------------------
jags.fit.mcmc %>%
  spread_draws(`beta_c[1,1]`, `beta_c[2,1]`, `beta_c[3,1]`, `beta_c[4,1]`, `beta_c[5,1]`,
               `beta_c[1,2]`, `beta_c[2,2]`, `beta_c[3,2]`, `beta_c[4,2]`, `beta_c[5,2]`,
               `beta_c[1,3]`, `beta_c[2,3]`, `beta_c[3,3]`, `beta_c[4,3]`, `beta_c[5,3]`,
               `beta_c[1,4]`, `beta_c[2,4]`, `beta_c[3,4]`, `beta_c[4,4]`, `beta_c[5,4]`) %>%
  pivot_longer(cols = c(`beta_c[1,1]`, `beta_c[2,1]`, `beta_c[3,1]`, `beta_c[4,1]`, `beta_c[5,1]`,
                        `beta_c[1,2]`, `beta_c[2,2]`, `beta_c[3,2]`, `beta_c[4,2]`, `beta_c[5,2]`,
                        `beta_c[1,3]`, `beta_c[2,3]`, `beta_c[3,3]`, `beta_c[4,3]`, `beta_c[5,3]`,
                        `beta_c[1,4]`, `beta_c[2,4]`, `beta_c[3,4]`, `beta_c[4,4]`, `beta_c[5,4]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  facet_wrap(~Parameter, ncol=4)
ggsave(filename = "betac_plot_all_bspline.pdf", path = vispath, height=12, width=15) 

jags.fit.mcmc %>%
  spread_draws(`a0[1,1]`, `a0[2,1]`, `a0[3,1]`, `a0[4,1]`, `a0[5,1]`, `a0[6,1]`,
               `a0[1,2]`, `a0[2,2]`, `a0[3,2]`, `a0[4,2]`, `a0[5,2]`, `a0[6,2]`,
               `a0[1,3]`, `a0[2,3]`, `a0[3,3]`, `a0[4,3]`, `a0[5,3]`, `a0[6,3]`,
               `a0[1,4]`, `a0[2,4]`, `a0[3,4]`, `a0[4,4]`, `a0[5,4]`, `a0[6,4]`,
               `a0[1,5]`, `a0[2,5]`, `a0[3,5]`, `a0[4,5]`, `a0[5,5]`, `a0[6,5]`
  ) %>%
  pivot_longer(cols = c(`a0[1,1]`, `a0[2,1]`, `a0[3,1]`, `a0[4,1]`, `a0[5,1]`, `a0[6,1]`,
                        `a0[1,2]`, `a0[2,2]`, `a0[3,2]`, `a0[4,2]`, `a0[5,2]`, `a0[6,2]`,
                        `a0[1,3]`, `a0[2,3]`, `a0[3,3]`, `a0[4,3]`, `a0[5,3]`, `a0[6,3]`,
                        `a0[1,4]`, `a0[2,4]`, `a0[3,4]`, `a0[4,4]`, `a0[5,4]`, `a0[6,4]`,
                        `a0[1,5]`, `a0[2,5]`, `a0[3,5]`, `a0[4,5]`, `a0[5,5]`, `a0[6,5]`
                        ), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  facet_wrap(~Parameter, ncol=5)
ggsave(filename = "a0_plot_all_bspline.pdf", path = vispath, height=12, width=15) 

# Posterior predictive check
y_sim <- fit_sim$BUGSoutput$sims.list$Y_sim
y_sim_mean <- apply(y_sim, 2, mean)
y_obs <- as.vector(unlist(logit.data[,c("logit.Public")]))

pdf(paste0(vispath, "posterior_predictive_observed_all_bspline.pdf")) 
plot(y_sim_mean, y_obs)
abline(a=0, b=1, col='red')
dev.off()

pdf(paste0(vispath,"posterior_predictive_check_all_bspline.pdf"))
# plot the observed data hist
hist(y_obs, freq = FALSE, breaks = 50) # observed data 
n_samples = nrow(y_sim) 
for (r in 1:n_samples){ 
  lines(density(y_sim[r,]), col = rgb(135, 206, 235, max= 255, alpha=25)) 
}
dev.off()

# Get P estimates --------------------------------------------------------------
# P1 <- readRDS('results/JAGS/JAGS_mod_Bspline_P_samps_chain1.RDS')
# P2 <- readRDS('results/JAGS/JAGS_mod_Bspline_P_samps_chain2.RDS')
# P3 <- readRDS('results/JAGS/JAGS_mod_Bspline_P_samps_chain3.RDS')
# 

# P_samps.mean <- apply(readRDS('results/JAGS/P_samps/Bspline/P_array_bspline_TG.RDS'), c(2,3,4,5), mean)
# P_samps.mean <- plyr::adply(P_samps.mean, .margins=c(1,2,4))
# colnames(P_samps.mean) <- c('index_subnat', 'index_method', 'index_year', 'Public', 'Private') 
# P_samps.mean <- P_samps.mean %>% 
#   mutate(across(everything(), as.numeric)) %>%
#   left_join(method_index_table) %>%
#   pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')
# 
# P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
# P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
# dim(P_samps.quantile)
# colnames(P_samps.quantile) <- c('index_subnat', 'index_method', 'index_sector', 'index_year', 'lower_95', 'upper_95') 
# P_samps.quantile <- P_samps.quantile %>% 
#   mutate(across(everything(), as.numeric)) %>%
#   left_join(method_index_table) %>%
#   left_join(sector_index_table) %>%
#   left_join(index_country_subnat_tbl) %>%
#   left_join(year_index_table)
# 
# P_samps_df <- left_join(P_samps.mean, P_samps.quantile)
# # saveRDS(P_samps_df, file='results/JAGS/P_samps/JAGS_bspline_4country_Psamps_df.RDS')

P_samps_df <- readRDS('results/JAGS/P_samps/Bspline/JAGS_mod_Bspline_TG_P_samps_df.RDS')

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
    theme(title = element_text(size=20), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
    theme(legend.position = "bottom", legend.title = element_text(size = 20), legend.text = element_text(size = 20))+
    facet_wrap(~interaction(Region, Method), ncol=5)
  ggsave(filename = paste0("p_plot_",c,".pdf"), path = vispath, height=12, width=15) 
  
}

# alpha_samps <- model_samps$a0
# dim(alpha_samps)
# alpha_mean <- t(apply(alpha_samps, c(2,3), mean))
# colnames(alpha_mean) <- n_subnat
# 
# alpha_mean <- as_tibble(alpha_mean) %>%
#   mutate(index_method = 1:length(n_method)) %>%
#   pivot_longer(cols=all_of(n_subnat), names_to = 'Region', values_to = 'alpha') %>%
#   mutate(invlogit.alpha = exp(alpha)/(1+exp(alpha))) %>%
#   left_join(method_index_table) %>%
#   left_join(index_country_subnat_tbl)
# 
# # Plot means vs observed values
# for(c in n_country) {
#   ggplot() +
#     geom_point(data = P_df %>% filter(Country==c), aes(x=average_year, y=Observed, colour=Sector, pch=Sector)) +
#     geom_errorbar(data = P_df %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax=upper_95, colour=Sector)) +
#     geom_line(data = P_samps_df %>% filter(Country==c), aes(x=average_year, y=Mean, colour=Sector, lty=Sector)) +
#     geom_ribbon(data = P_samps_df %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
#     geom_hline(data=alpha_mean %>% filter(Country==c), aes(yintercept = invlogit.alpha)) +
#     theme(title = element_text(size=20), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
#     theme(legend.position = "bottom", legend.title = element_text(size = 20), legend.text = element_text(size = 20))+
#     facet_wrap(~interaction(Method, Region), ncol=5)
#   ggsave(filename = paste0("p_plot_",c,"_withAlpha.pdf"), path = vispath, height=12, width=15) 
# }
