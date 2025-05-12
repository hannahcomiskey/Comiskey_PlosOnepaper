library(rstan)
library(tidyverse)
library(tidybayes)
library(bayesplot)

source('code/load_functions.R')
source('code/stan_utility.R')

# Fit model -----------------------------------------
vispath = "visualisations/simulated_data/bspline/STAN/10years/"
  
fit <- readRDS('results/STAN_model_basis_expand_KenyaSim_10_cauchy.RDS')

code <- get_stancode(fit)
cat(code)

check_all_diagnostics(fit)

# pairs(fit, pars = "a0")

check_divergences(fit)

# Get parameter estimates 
model_samps <- rstan::extract(fit)

fit_mcmc <- coda::as.mcmc(fit)

# Traceplots
traceplot(fit, pars = c("a[1,5,1]", "a[2,5,3]", "a[3,5,9]", "a[4,5,6]", "a[5,5,4]", "a[6,5,4]"), inc_warmup = FALSE, nrow=6)
traceplot(fit, pars = c("a0"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("tau"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_y"), inc_warmup = FALSE, nrow = 5)


# Check divergences 
c_dark <- c("#8F272780")
green <- c("#00FF0080")

partition <- partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mar = c(4, 4, 0.5, 0.5))
plot(nondiv_params$`a0[1,1]`, nondiv_params$`sigma_a0[1]`,
     col=c_dark, pch=16, cex=0.8)
points(div_params$`a0[1,1]`, div_params$`sigma_a0[1]`,
       col='blue', pch=16, cex=0.8)

# Checking for funnel in beta parameters
par(mfrow=c(5, 1))
for (k in 1:length(n_method)) {
  name <- paste("a0[1,", k, "]", sep="")
  # name2 <- paste("sigma_a0[", k, "]", sep="")
  plot(nondiv_params[name][[1]], #log(nondiv_params[name2][[1]]),
       col="#8F272780", pch=16, cex=0.8,
       xlab=name) #, ylab="log(tau)")
  points(div_params[name][[1]], #log(div_params[name2][[1]]),
         col="blue", pch=16, cex=0.8)
}
par(mfrow=c(1, 1))


# Rhats 
rhats <- rhat(fit)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats) + yaxis_text(hjust = 1)


# N-eff
ratios <- bayesplot::neff_ratio(fit)
mcmc_neff(ratios)

mcmc_nuts_divergence(nuts_params(fit), log_posterior(fit))

np_cp <- nuts_params(fit)
mcmc_nuts_energy(np_cp)

fit_summary <- summary(fit)

# Review spline coefficients
beta.k_samps <- model_samps$a
beta_k.mean <- apply(beta.k_samps, c(2,3,4), mean)
dim(beta_k.mean)
sum(beta_k.mean[1,1,])

# Index tables ---------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))

# # Check priors vs posteriors samples -----------------
# beta_c_new <- t(as.matrix(beta_c_new))
# colnames(beta_c_new) <- 1:5
# beta_c_df <- as_tibble(beta_c_new) %>%
#   pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
#   mutate(index_country = rep(1, each=1),
#          Parameter = paste0("beta_c[", rep(1:5, 1),"]"))
# 
# fit %>%
#   spread_draws(`beta_c[1]`, `beta_c[2]`, `beta_c[3]`, `beta_c[4]`, `beta_c[5]`) %>%
#   pivot_longer(cols = c(`beta_c[1]`, `beta_c[2]`, `beta_c[3]`, `beta_c[4]`, `beta_c[5]`), names_to='Parameter', values_to = 'sample') %>%
#   ggplot() +
#   #stat_halfeye(aes(y = Parameter, x = sample)) +
#   geom_density(aes(x = sample), fill='grey', alpha=0.8) +
#   geom_vline(data=beta_c_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
#   facet_wrap(~Parameter)
# ggsave(filename = paste0(vispath, 'beta_c_density_KenyaSim50.pdf'))

# fit %>%
#   spread_draws(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`) %>%
#   pivot_longer(cols = c(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`), names_to='Parameter', values_to = 'sample') %>%
#   ggplot() +
#   #stat_halfeye(aes(y = Parameter, x = sample)) +
#   geom_density(aes(x = sample), fill='grey', alpha=0.8) +
#   geom_vline(data=beta_c_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
#   facet_wrap(~Parameter)
# ggsave(filename = paste0(vispath, 'beta_c_density_KenyaSim50.pdf'))


# Get P estimates 

P_samps <- model_samps$P
dim(P_samps)
P_samps.mean <- apply(P_samps, c(2,3,4,5), mean)
P_samps.mean <- plyr::adply(P_samps.mean, .margins=c(1,2,4))
colnames(P_samps.mean) <- c('index_subnat', 'index_method',  'index_year', 'Public', 'Private') 
P_samps.mean <- P_samps.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')

P_samps <- model_samps$P
dim(P_samps)
P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
colnames(P_samps.quantile) <- c('index_subnat', 'index_method', 'index_sector', 'index_year', 'lower_95', 'upper_95') 
P_samps.quantile <- P_samps.quantile %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  left_join(sector_index_table) 

P_samps_df <- left_join(P_samps.mean, P_samps.quantile)

# Get observed data 
P_df<- P_sim_df_sample %>%
  mutate(across(everything(), as.numeric)) %>%
  left_join(year_index_table) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

# Plot means vs observed values
ggplot() +
  geom_point(data = P_df, aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_samps_df, aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df, aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  facet_wrap(~interaction(Method, index_subnat), ncol=5)
ggsave(filename = "p_plot_Kenyasim_10years.pdf", path = vispath, height=12, width=15) 


