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
vispath = 'visualisations/MVN_delta/'

# fit_sim model -----------------------------------------
fit_sim <- readRDS('results/JAGS/JAGS_model_MVN_NCP_kstar_SE_sum0_dsigma_parallel_noindia.RDS')

fit_sim$model

plot(fit_sim)

model_samps <- fit_sim$BUGSoutput$sims.list

View(fit_sim$BUGSoutput$summary)

jags.fit.mcmc <- as.mcmc(fit_sim)

# Plot ACF ---------------------------------------------------------------------
mcmc_acf(jags.fit.mcmc, pars='sigma_delta')
mcmc_acf(jags.fit.mcmc, pars=c('sigma_alpha_cms[1]'))
mcmc_acf(jags.fit.mcmc, pars=c('sigma_alpha_pms[1]'))

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))

# Plot traceplots -------------------------------------------------------------- 
mcmc_trace(jags.fit.mcmc, pars= c('alpha_cms[1,1]', 'alpha_cms[2,1]', 'alpha_cms[3,1]', 'alpha_cms[4,1]', 'alpha_cms[5,1]'))

mcmc_trace(jags.fit.mcmc, pars= c('err_alpha_cms[1,1]', 'err_alpha_cms[2,1]', 'err_alpha_cms[3,1]', 'err_alpha_cms[4,1]', 'err_alpha_cms[5,1]'))
mcmc_trace(jags.fit.mcmc, pars= c('err_alpha_cms[1,1]', 'err_alpha_cms[1,2]', 'err_alpha_cms[1,3]', 'err_alpha_cms[1,4]', 'err_alpha_cms[1,5]'))

mcmc_trace(jags.fit.mcmc, pars= c('err_alpha_pms[1,1,1]', 'err_alpha_pms[2,1,2]', 'err_alpha_pms[3,1,3]', 'err_alpha_pms[4,1,4]', 'err_alpha_pms[5,1,5]'))
mcmc_trace(jags.fit.mcmc, pars= c('err_alpha_pms[1,1,1]', 'err_alpha_pms[2,2,2]', 'err_alpha_pms[3,3,3]', 'err_alpha_pms[4,4,4]', 'err_alpha_pms[5,5,5]'))

mcmc_trace(jags.fit.mcmc, pars= c('sigma_delta'))

# Check priors vs posteriors samples -------------------------------------------
jags.fit.mcmc %>%
  spread_draws(`alpha_cms[1,1]`, `alpha_cms[2,1]`, `alpha_cms[3,1]`, `alpha_cms[4,1]`, `alpha_cms[5,1]`,
               `alpha_cms[1,2]`, `alpha_cms[2,2]`, `alpha_cms[3,2]`, `alpha_cms[4,2]`, `alpha_cms[5,2]`,
               `alpha_cms[1,3]`, `alpha_cms[2,3]`, `alpha_cms[3,3]`, `alpha_cms[4,3]`, `alpha_cms[5,3]`,
               `alpha_cms[1,4]`, `alpha_cms[2,4]`, `alpha_cms[3,4]`, `alpha_cms[4,4]`, `alpha_cms[5,4]`) %>%
  pivot_longer(cols = c(`alpha_cms[1,1]`, `alpha_cms[2,1]`, `alpha_cms[3,1]`, `alpha_cms[4,1]`, `alpha_cms[5,1]`,
                        `alpha_cms[1,2]`, `alpha_cms[2,2]`, `alpha_cms[3,2]`, `alpha_cms[4,2]`, `alpha_cms[5,2]`,
                        `alpha_cms[1,3]`, `alpha_cms[2,3]`, `alpha_cms[3,3]`, `alpha_cms[4,3]`, `alpha_cms[5,3]`,
                        `alpha_cms[1,4]`, `alpha_cms[2,4]`, `alpha_cms[3,4]`, `alpha_cms[4,4]`, `alpha_cms[5,4]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  facet_wrap(~Parameter, ncol=4)
ggsave(filename = "betac_plot_.pdf", path = paste0(vispath,'/parameter_plots/'), height=12, width=15) 


# Posterior predictive check
y_sim <- fit_sim$BUGSoutput$sims.list$y.sim
y_sim_mean <- apply(y_sim, 2, mean)
y_obs <- as.vector(unlist(logit.data[,c("logit.Public")]))

pdf(paste0(vispath, "posterior_predictive_observed_MVN.pdf")) 
plot(y_sim_mean, y_obs)
abline(a=0, b=1, col='red')
dev.off()

pdf(paste0(vispath,"posterior_predictive_check_MVN.pdf"))
# plot the observed data hist
hist(y_obs, freq = FALSE, breaks = 50) # observed data 
n_samples = nrow(y_sim) 
for (r in 1:n_samples){ 
  lines(density(y_sim[r,]), col = rgb(135, 206, 235, max= 255, alpha=25)) 
}
dev.off()

# Get P estimates --------------------------------------------------------------

# P_samps <- model_samps$P
# dim(P_samps)
# P_samps.mean <- apply(P_samps, c(2,3,4,5), mean)
# dim(P_samps.mean)
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
# colnames(P_samps.quantile) <- c( 'index_sector', 'index_method', 'index_subnat', 'index_year', 'lower_95', 'upper_95') 
# P_samps.quantile <- P_samps.quantile %>% 
#   mutate(across(everything(), as.numeric)) %>%
#   left_join(method_index_table) %>%
#   left_join(sector_index_table) %>%
#   left_join(index_country_subnat_tbl) %>%
#   left_join(year_index_table)
# 
# P_samps_df <- left_join(P_samps.mean, P_samps.quantile)
# saveRDS(P_samps_df, file='results/JAGS/P_samps/JAGS_model_all_MVN_Psamps_df.RDS')

P_samps_df <- readRDS('results/JAGS/P_samps/MVN_delta_chains/JAGS_MVN_delta_Psamps_df.RDS')

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
  ggsave(filename = paste0("p_plot_",c,".pdf"), path = paste0(vispath, '/country_plots'), height=12, width=15) 
  
}

alpha_samps <- model_samps$alpha_pms
dim(alpha_mean)
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
    theme(title = element_text(size=20), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=8), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
    theme(legend.position = "bottom", legend.title = element_text(size = 20), legend.text = element_text(size = 20))+
    facet_wrap(~interaction(Method, Region), ncol=5)
  ggsave(filename = paste0("p_plot_",c,"_withAlpha.pdf"), path = paste0(vispath, '/country_plots'), height=12, width=15) 
}
