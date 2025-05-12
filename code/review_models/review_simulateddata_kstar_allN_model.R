library(rstan)
library(tidyverse)
library(tidybayes)
source('code/stan_utility.R')
source('code/load_functions.R')

vispath= 'visualisations/simulated_data/kstar/all_N/'

# fit_sim model -----------------------------------------
fit_sim <- readRDS('results/STAN_model_allN_NCP_kstar_sum0_KenyaSim_10.RDS')

code_sim <- get_stancode(fit_sim)
cat(code_sim)

check_all_diagnostics(fit_sim)

check_divergences(fit_sim)

print(fit_sim)

# Get parameter estimates 
model_samps <- rstan::extract(fit_sim)

# Traceplots
traceplot(fit_sim, pars = c("alpha_pms[3,1]", "alpha_pms[3,2]", "alpha_pms[3,3]", "alpha_pms[3,4]", "alpha_pms[3,5]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit_sim, pars = c("delta_k[1,1,1]",  "delta_k[2,3,5]",  "delta_k[3,2,2]", "delta_k[4,5,7]", "delta_k[5,4,6]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit_sim, pars = c("beta_c"), inc_warmup = FALSE, nrow = 5)
traceplot(fit_sim, pars = c("beta_k[1,1,1]",  "beta_k[2,3,5]",  "beta_k[3,2,2]", "beta_k[4,5,7]", "beta_k[5,4,6]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit_sim, pars = c("sigma_delta"), inc_warmup = FALSE, nrow = 5)
traceplot(fit_sim, pars = c("sigma_beta"), inc_warmup = FALSE, nrow = 5)
traceplot(fit_sim, pars = c("sigma_alpha"), inc_warmup = FALSE, nrow = 5)


# Check divergences 
c_dark <- c("#8F272780")
green <- c("#00FF0080")

partition <- partition_div(fit_sim)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mar = c(4, 4, 0.5, 0.5))
plot(nondiv_params$`beta_c[1,2]`, nondiv_params$`sigma_beta[1]`,
     col=c_dark, pch=16, cex=0.8)
points(div_params$`beta_c[1,1]`, div_params$`sigma_beta[1,1]`,
       col='blue', pch=16, cex=0.8)

# Checking for funnel in beta parameters
par(mfrow=c(5, 1))
for (k in 1:length(n_method)) {
  name <- paste("beta_c[1,", k, "]", sep="")
  name2 <- "sigma_beta" #paste("sigma_beta[", k, "]", sep="")
  plot(nondiv_params[name][[1]], log(nondiv_params[name2][[1]]),
       col="#8F272780", pch=16, cex=0.8,
       xlab=name, ylab="log(tau)")
  points(div_params[name][[1]], log(div_params[name2][[1]]),
         col="blue", pch=16, cex=0.8)
}
par(mfrow=c(1, 1))


# Rhats 
rhats <- rhat(fit_sim)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats) + yaxis_text(hjust = 1)

# N-eff
ratios <- bayesplot::neff_ratio(fit_sim)
mcmc_neff(ratios)

mcmc_nuts_divergence(nuts_params(fit_sim), log_posterior(fit_sim))

np_cp <- nuts_params(fit_sim)
mcmc_nuts_energy(np_cp)

fit_sim_summary <- summary(fit_sim)

# Check delta.k sums to 0 
beta.k_samps <- model_samps$beta_k
dim(beta.k_samps)
sum(beta.k_samps[1,1,1,]) + model_samps$alpha_pms[1,1,1]

# Review spline coefficients
beta_k.mean <- apply(beta.k_samps, c(2,3,4), mean)

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(Year = all_years, index_year = 1:length(all_years))

# Check priors vs posteriors samples -------------------------------------------
beta_c_new <- t(as.matrix(beta_c_new))
colnames(beta_c_new) <- 1:5
beta_c_df <- as_tibble(beta_c_new) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  mutate(index_country = rep(1, each=1),
         Parameter = paste0("beta_c[",rep(1, each=5),",",rep(1:5, 1),"]"))

fit_sim %>%
  spread_draws(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`) %>%
               # `beta_c[2,1]`, `beta_c[2,2]`, `beta_c[2,3]`, `beta_c[2,4]`, `beta_c[2,5]`,
               # `beta_c[3,1]`, `beta_c[3,2]`, `beta_c[3,3]`, `beta_c[3,4]`, `beta_c[3,5]`,
               # `beta_c[4,1]`, `beta_c[4,2]`, `beta_c[4,3]`, `beta_c[4,4]`, `beta_c[4,5]`
  pivot_longer(cols = c(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  #stat_halfeye(aes(y = Parameter, x = sample)) +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  geom_vline(data=beta_c_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = "betac_plot_Kenyasim_10years.pdf", path = vispath, height=12, width=15) 


alpha_sim <- t(alpha_sim)
colnames(alpha_sim) <- 1:5
alpha_df <- as_tibble(alpha_sim) %>%
  mutate(index_subnat = 1:P) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  rowwise() %>%
  mutate(Parameter = paste0("alpha_pms[", index_subnat,",",index_method,"]"))
  # mutate(Parameter = paste0("alpha_pms[",index_subnat,",",index_method,"]"))

fit_sim %>%
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
  geom_vline(data=alpha_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = "alpha_plot_Kenyasim_10years.pdf", path = vispath, height=12, width=15)

# fit_sim %>%
#   spread_draws(`alpha_pms[1,1]`, `alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[6,1]`,
#                `alpha_pms[1,2]`, `alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[6,2]`,
#                `alpha_pms[1,3]`, `alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[6,3]`,
#                `alpha_pms[1,4]`, `alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[6,4]`,
#                `alpha_pms[1,5]`, `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[6,5]`
#   ) %>%
#   pivot_longer(cols = c(`alpha_pms[1,1]`, `alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[6,1]`,
#                         `alpha_pms[1,2]`, `alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[6,2]`,
#                         `alpha_pms[1,3]`, `alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[6,3]`,
#                         `alpha_pms[1,4]`, `alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[6,4]`,
#                         `alpha_pms[1,5]`, `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[6,5]`), names_to='Parameter', values_to = 'sample') %>%
#   ggplot() +
#   geom_density(aes(x = sample), fill='grey', alpha=0.8) +
#   geom_vline(data=alpha_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
#   facet_wrap(~Parameter)
# ggsave(filename = 'visualisations/simulated_data/kstar/all_N/ncp_alpha_pm_density_KenyaSim_10years.pdf')


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
  left_join(sector_index_table) 

P_samps_df <- left_join(P_samps.mean, P_samps.quantile)
# saveRDS(P_samps_df, file='results/STAN_model_Kenya_allN_NCP_kstar_10years_Psamps_df.RDS')

# Get observed data 
P_df<- P_sim_df_sample %>%
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

#%>% filter(index_year<30)
# Plot means vs observed values
ggplot() +
  geom_point(data = P_df, aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_samps_df, aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df , aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  facet_wrap(~interaction(Method, index_subnat), ncol=5)
ggsave(filename = "p_plot_Kenyasim_10years.pdf", path = vispath, height=12, width=15) 


alpha_samps <- model_samps$alpha_pms
dim(alpha_samps)
alpha_mean <- t(apply(alpha_samps, c(2,3), mean))
colnames(alpha_mean) <- 1:P

alpha_mean <- as_tibble(alpha_mean) %>%
  mutate(index_method = 1:M) %>%
  pivot_longer(cols=`1`:`6`, names_to = 'index_subnat', values_to = 'alpha') %>%
  mutate(invlogit.alpha = exp(alpha)/(1+exp(alpha))) %>%
  left_join(method_index_table)

ggplot() +
  geom_point(data = P_df, aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_samps_df, aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df , aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  geom_hline(data=alpha_mean, aes(yintercept = invlogit.alpha), lty='dashed') +
  facet_wrap(~interaction(Method, index_subnat), ncol = 5)
ggsave(filename = "p_plot_Kenyasim_10years_withAlpha.pdf", path = vispath, height=12, width=15) 
