library(rstan)
library(tidyverse)
library(tidybayes)
library(bayesplot)
source('code/load_functions.R')
source('code/stan_utility.R')

# Fit model -----------------------------------------
fit <- readRDS('results/STAN_model_allN_KenyaSim_50.RDS')

code <- get_stancode(fit)
cat(code)

traceplot(fit, pars = c("alpha_pms[3,1]", "alpha_pms[3,2]", "alpha_pms[3,3]", "alpha_pms[3,4]", "alpha_pms[3,5]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("delta_k[1,1,1]",  "delta_k[2,3,5]",  "delta_k[3,2,2]", "delta_k[4,5,7]", "delta_k[5,4,6]"), inc_warmup = FALSE, nrow = 6)
# traceplot(fit, pars = c("sigmabeta_tau"), inc_warmup = FALSE, nrow = 6)
# traceplot(fit, pars = c("sigmabeta_Omega"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("sigma_alpha"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_delta"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("beta_c"), inc_warmup = FALSE, nrow = 5)

pairs(fit, pars = "sigma_y")

check_divergences(fit)

check_all_diagnostics(fit)

# Checking for funnel in beta parameters
partition <- partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]
par(mfrow=c(5, 1))
for (k in 1:length(n_method)) {
  name <- paste("beta_c[1,", k, "]", sep="")
  name2 <- paste("sigma_beta[", k, "]", sep="")
  plot(nondiv_params[name][[1]], log(nondiv_params[name2][[1]]),
       col="#8F272780", pch=16, cex=0.8,
       xlab=name, xlim=c(-2.5, 2.5), ylab="log(tau)")
  points(div_params[name][[1]], log(div_params[name2][[1]]),
         col="blue", pch=16, cex=0.8)
}
par(mfrow=c(1, 1))

# Rhats 
rhats <- rhat(fit)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats) + yaxis_text(hjust = 1)


# N-eff
ratios <- bayesplot::neff_ratio(fit)
print(ratios)
mcmc_neff(ratios)

mcmc_nuts_divergence(nuts_params(fit), log_posterior(fit))

np_cp <- nuts_params(fit)
mcmc_nuts_energy(np_cp)

fit_summary <- summary(fit)
View(print(fit_summary$summary))

# Get parameter estimates 

model_samps <- rstan::extract(fit)

# Check delta.k sums to 0 
delta.k_samps <- model_samps$delta_k
dim(delta.k_samps)
sum(delta.k_samps[1,1,1,1:12])

# Review spline coefficients
delta_k.mean <- apply(delta.k_samps, c(2,3,4), mean)

# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- matchcountry
n_all_years <- length(all_years)

# Check out variance-covariance matrices 
View(get_posterior_mean(fit))

# Compare posterior to inputs 

get_variables(fit)

beta_c_new <- t(as.matrix(beta_c_new))
colnames(beta_c_new) <- 1:5
beta_c_df <- as_tibble(beta_c_new) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  mutate(index_country = rep(1, each=1),
         Parameter = paste0("beta_c[1,", rep(1:5, 1),"]"))

fit %>%
  spread_draws(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`) %>%
  pivot_longer(cols = c(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  #stat_halfeye(aes(y = Parameter, x = sample)) +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  geom_vline(data=beta_c_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = 'visualisations/simulated_data/kstar/all_N/beta_c_density.pdf')


#HERE
alpha_sim <- t(alpha_sim)
colnames(alpha_sim) <- 1:5
alpha_df <- as_tibble(alpha_sim) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  mutate(index_subnat = rep(1:P, each=5),
         Parameter = paste0("alpha_pms[",rep(1:M, M),",",rep(1:P, each=M),"]"))

fit %>%
  spread_draws(`alpha_pms[1,1]`, `alpha_pms[1,2]`, `alpha_pms[1,3]`, `alpha_pms[1,4]`, `alpha_pms[1,5]`, `alpha_pms[1,6]`,
               `alpha_pms[2,1]`, `alpha_pms[2,2]`, `alpha_pms[2,3]`, `alpha_pms[2,4]`, `alpha_pms[2,5]`, `alpha_pms[2,6]`,
               `alpha_pms[3,1]`, `alpha_pms[3,2]`, `alpha_pms[3,3]`, `alpha_pms[3,4]`, `alpha_pms[3,5]`, `alpha_pms[3,6]`,
               `alpha_pms[4,1]`, `alpha_pms[4,2]`, `alpha_pms[4,3]`, `alpha_pms[4,4]`, `alpha_pms[4,5]`, `alpha_pms[4,6]`,
               `alpha_pms[5,1]`, `alpha_pms[5,2]`, `alpha_pms[5,3]`, `alpha_pms[5,4]`, `alpha_pms[5,5]`, `alpha_pms[5,6]`) %>%
  pivot_longer(cols = c(`alpha_pms[1,1]`, `alpha_pms[1,2]`, `alpha_pms[1,3]`, `alpha_pms[1,4]`, `alpha_pms[1,5]`, `alpha_pms[1,6]`,
                        `alpha_pms[2,1]`, `alpha_pms[2,2]`, `alpha_pms[2,3]`, `alpha_pms[2,4]`, `alpha_pms[2,5]`, `alpha_pms[2,6]`,
                        `alpha_pms[3,1]`, `alpha_pms[3,2]`, `alpha_pms[3,3]`, `alpha_pms[3,4]`, `alpha_pms[3,5]`, `alpha_pms[3,6]`,
                        `alpha_pms[4,1]`, `alpha_pms[4,2]`, `alpha_pms[4,3]`, `alpha_pms[4,4]`, `alpha_pms[4,5]`, `alpha_pms[4,6]`,
                        `alpha_pms[5,1]`, `alpha_pms[5,2]`, `alpha_pms[5,3]`, `alpha_pms[5,4]`, `alpha_pms[5,5]`, `alpha_pms[5,6]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  geom_vline(data=alpha_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = 'visualisations/simulated_data/kstar/all_N/alpha_pm_density.pdf')


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

sector_index_table <- tibble(Sector = c('Public', 'Private'), index_sector = 1:2)
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

# Get alpha intercepts 
alpha_samps <- model_samps$alpha_pms
dim(alpha_samps)
alpha_samps.mean <- as_tibble(apply(alpha_samps, c(2,3), mean)) 
colnames(alpha_samps.mean) <- c(1:P)
alpha_samps.mean <- alpha_samps.mean %>%
  mutate(index_method = rep(1:5)) %>%
  pivot_longer(cols=`1`:`6`, names_to = 'index_subnat', values_to = 'alpha_mean')
alpha_samps.mean <- alpha_samps.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  mutate(invlogit.alpha = exp(alpha_mean)/(1+exp(alpha_mean)))

apply(model_samps$sigma_alpha, 2, mean)

# Get observed data 
P_df<- P_sim_df_sample %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

# Plot means vs observed values
ggplot() +
  geom_point(data = P_df, aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_samps_df, aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df, aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  facet_wrap(~interaction(Method, index_subnat), ncol=5)
ggsave(filename = 'visualisations/simulated_data/kstar/all_N/plotting_p.pdf')



