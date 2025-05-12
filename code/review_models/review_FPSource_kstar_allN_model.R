library(rstan)
library(tidyverse)
library(tidybayes)
source('code/stan_utility.R')
source('code/load_functions.R')


# Source code --------------------------------------
source("code/load_functions.R")
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

vispath = 'visualisations/country_data/allN/STAN/'

# Get logit of parameters and variance --------------------
mydata <- FP_source_data_wide[,c("Public", "Public.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))

# Fit model -----------------------------------------
fit <- readRDS(file='results/STAN/STAN_model_4country_allN_NCP_kstar.RDS')

code_sim <- get_stancode(fit)
cat(code_sim)

check_energy(fit)
check_divergences(fit)
check_n_eff(fit)
check_rhat(fit)
check_treedepth(fit)

stan_ac(fit, pars = c('sigma_alpha', 'sigma_delta', 'sigma_beta'))
stan_ac(fit, pars = 'beta_c')

stan_plot(fit,  pars = c('sigma_alpha', 'sigma_delta', 'sigma_beta'))

stan_hist(fit,  pars = c('sigma_alpha', 'sigma_delta', 'sigma_beta'))

stan_dens(fit,  pars = c('sigma_alpha', 'sigma_delta', 'sigma_beta'))

print(fit)

# Get parameter estimates 
model_samps <- rstan::extract(fit)

# Traceplots
traceplot(fit, pars = c("alpha_pms[3,1]", "alpha_pms[3,2]", "alpha_pms[3,3]", "alpha_pms[3,4]", "alpha_pms[3,5]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("delta_k[1,1,1]",  "delta_k[2,3,5]",  "delta_k[3,2,2]", "delta_k[4,5,7]", "delta_k[5,4,6]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("beta_c"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("beta_k[1,1,1]",  "beta_k[2,3,5]",  "beta_k[3,2,2]", "beta_k[4,5,7]", "beta_k[5,4,6]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("sigma_delta"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_beta"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_alpha"), inc_warmup = FALSE, nrow = 5)

code <- get_stancode(fit)
cat(code)

check_divergences(fit)

# Rhats 
rhats <-  bayesplot::rhat(fit)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats) + yaxis_text(hjust = 1)

# N-eff
ratios <- bayesplot::neff_ratio(fit)
mcmc_neff(ratios)

mcmc_nuts_divergence(nuts_params(fit), log_posterior(fit))

np_cp <- nuts_params(fit)
mcmc_nuts_energy(np_cp)

# Get parameter estimates 
model_samps <- rstan::extract(fit)

fit_summary <- summary(fit)

# Check delta.k sums to 0 
beta.k_samps <- model_samps$beta_k
dim(beta.k_samps)
sum(beta.k_samps[1,1,1,]) + model_samps$alpha_pms[1,1,1]

# Review spline coefficients
beta_k.mean <- apply(beta.k_samps, c(2,3,4), mean)

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))

# Check priors vs posteriors samples -------------------------------------------
fit %>%
  spread_draws(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`) %>%
  pivot_longer(cols = c(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  #stat_halfeye(aes(y = Parameter, x = sample)) +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  facet_wrap(~Parameter)
ggsave(filename = "betac_plot_4country.pdf", path = vispath, height=12, width=15) 


fit %>%
  spread_draws(`alpha_pms[1,1]`, `alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[6,1]`,
               `alpha_pms[1,2]`, `alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[6,2]`,
               `alpha_pms[1,3]`, `alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[6,3]`,
               `alpha_pms[1,4]`, `alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[6,4]`,
               `alpha_pms[1,5]`, `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[6,5]`
  ) %>%
  pivot_longer(cols = c(`alpha_pms[1,1]`, `alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[6,1]`,
                        `alpha_pms[1,2]`, `alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[6,2]`,
                        `alpha_pms[1,3]`, `alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[6,3]`,
                        `alpha_pms[1,4]`, `alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[6,4]`,
                        `alpha_pms[1,5]`, `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[6,5]`
  ), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  facet_wrap(~Parameter)
ggsave(filename = "alpha_plot_4country.pdf", path = vispath, height=12, width=15)


# Get P estimates 

P_samps <- model_samps$P
dim(P_samps)
P_samps.mean <- apply(P_samps, c(2,3,4,5), mean)
P_samps.mean <- plyr::adply(P_samps.mean, .margins=c(1,2,4))
colnames(P_samps.mean) <- c('index_method', 'index_subnat', 'index_year', 'Public', 'Private') 
P_samps.mean <- P_samps.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')

P_samps <- model_samps$P
dim(P_samps)
P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
colnames(P_samps.quantile) <- c('index_method', 'index_subnat','index_sector', 'index_year', 'lower_95', 'upper_95') 
P_samps.quantile <- P_samps.quantile %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  left_join(sector_index_table) %>%
  left_join(index_country_subnat_tbl) %>%
  left_join(year_index_table)

P_samps_df <- left_join(P_samps.mean, P_samps.quantile)
saveRDS(P_samps_df, file='results/STAN_model_4country_allN_NCP_Psamps_df.RDS')

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
